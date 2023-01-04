import errno
import os
import tarfile
from subprocess import Popen, PIPE

import csverve.api as csverve
import pandas as pd
import yaml
from remixt.analysis.haplotype import calculate_haplotypes
from remixt.analysis.haplotype import read_phasing_samples


def makedirs(directory, isfile=False):
    if isfile:
        directory = os.path.dirname(directory)
        if not directory:
            return

    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def untar(input_tar, outdir):
    makedirs(outdir)
    with tarfile.open(input_tar) as tar:
        tar.extractall(path=outdir)


def load_tar_meta(thousand_genomes_dir):
    tg_tar_yaml = os.path.join(thousand_genomes_dir, 'meta.yaml')
    assert os.path.exists(tg_tar_yaml)

    with open(tg_tar_yaml, 'rt') as reader:
        data = yaml.safe_load(reader)

    return data


def run_cmd(cmd, output=None):
    cmd = [str(v) for v in cmd]

    print(' '.join(cmd))

    stdout = PIPE
    if output:
        stdout = open(output, "w")

    p = Popen(cmd, stdout=stdout, stderr=PIPE)

    cmdout, cmderr = p.communicate()
    retc = p.returncode

    if retc:
        raise Exception(
            "command failed. stderr:{}, stdout:{}".format(
                cmdout,
                cmderr))

    if output:
        stdout.close()


def run_shapeit_grch38(
        thousand_genomes_dir, chromosome, input_bcf, output_csv, temp_directory, shapeit_num_samples=100,
        shapeit_confidence_threshold=0.95
):
    meta = load_tar_meta(thousand_genomes_dir)

    assert '{chromosome}' in meta['grch38_1kg_bcf_filename_template']
    assert '{chromosome}' in meta['genetic_map_grch38_filename_template']

    grch38_genetic_map = meta['genetic_map_grch38_filename_template'].replace('{chromosome}', chromosome)
    grch38_genetic_map = os.path.join(thousand_genomes_dir, grch38_genetic_map)

    if chromosome == meta['grch38_1kg_phased_chromosome_x']:
        grch38_1kg_bcf = os.path.join(thousand_genomes_dir, meta['grch38_1kg_X_bcf_filename_template'])
    else:
        grch38_1kg_bcf = meta['grch38_1kg_bcf_filename_template'].replace('{chromosome}', chromosome)
        grch38_1kg_bcf = os.path.join(thousand_genomes_dir, grch38_1kg_bcf)

    bingraph_filename = os.path.join(temp_directory, 'phasing.bingraph')

    run_cmd(
        ['shapeit4',
         '--input', input_bcf,
         '--map', grch38_genetic_map,
         '--region', chromosome,
         '--reference', grch38_1kg_bcf,
         '--bingraph', bingraph_filename]
    )

    # Run shapeit to sample from phased haplotype graph
    sample_template = os.path.join(temp_directory, 'sampled.{0}.bcf')
    sample_filenames = []

    for s in range(shapeit_num_samples):
        sample_filename = sample_template.format(s)
        sample_filenames.append(sample_filename)

        run_cmd(
            ['bingraphsample',
             '--input', bingraph_filename,
             '--output', sample_filename,
             '--sample',
             '--seed', str(s)
             ]
        )

        run_cmd(
            ['bcftools', 'index', '-f', sample_filename]
        )

    haplotypes = calculate_haplotypes(
        read_phasing_samples(sample_filenames),
        changepoint_threshold=shapeit_confidence_threshold
    )

    haplotypes = pd.concat([
        haplotypes.rename(columns={'allele1': 'allele'})[['chromosome', 'position', 'allele', 'hap_label']].assign(
            allele_id=0),
        haplotypes.rename(columns={'allele2': 'allele'})[['chromosome', 'position', 'allele', 'hap_label']].assign(
            allele_id=1),
    ])

    # Translate from grch38 thousand genomes chr prefix
    if meta['chr_name_prefix'] == '':
        if not haplotypes['chromosome'].str.startswith('chr').all():
            raise ValueError('unexpected chromosome prefix')
        haplotypes['chromosome'] = haplotypes['chromosome'].str.slice(start=3)

    haplotypes = haplotypes[['chromosome', 'position', 'allele', 'hap_label', 'allele_id']]

    dtypes = {
        'chromosome': str,
        'position': int,
        'allele': int,
        'hap_label': int,
        'allele_id': int
    }

    csverve.write_dataframe_to_csv_and_yaml(haplotypes, output_csv, dtypes)


def infer_haps(
        input_bcf_file,
        output_csv,
        thousand_genomes_tar,
        chromosome,
        tempdir,
        shapeit_num_samples=100,
        shapeit_confidence_threshold=0.95
):
    thousand_genomes_dir = os.path.join(tempdir, 'thousand_genomes_impute_tar')

    untar(thousand_genomes_tar, thousand_genomes_dir)

    meta = load_tar_meta(thousand_genomes_dir)

    if meta['ensembl_genome_version'] == 'GRCh37':
        raise NotImplementedError
    elif meta['ensembl_genome_version'] == 'GRCh38':
        run_shapeit_grch38(
            thousand_genomes_dir,
            chromosome,
            input_bcf_file,
            output_csv,
            tempdir,
            shapeit_num_samples=shapeit_num_samples,
            shapeit_confidence_threshold=shapeit_confidence_threshold
        )
    else:
        raise Exception()
