import errno
import os
import tarfile
from subprocess import Popen, PIPE

import csverve.api as csverve
import numpy as np
import pandas as pd
import vcf
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

    chromosome = meta['grch38_1kg_phased_chromosome_x'] if chromosome.endswith('X') else chromosome

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


def write_null(haps_filename):
    with open(haps_filename, 'w') as haps_file:
        haps_file.write('chromosome\tposition\tallele\thap_label\tallele_id\n')


def read_bcf(bcf_file, chromosome):
    reader = vcf.Reader(filename=bcf_file)

    data = []
    for record in reader:
        chrom = record.CHROM
        pos = record.POS

        assert chrom == chromosome

        assert len(record.samples) == 1

        for sample in record.samples:
            gt = sample['GT']
            break

        if gt in ['0/0', '0|0']:
            AA = 1
            AB = 0
            BB = 0
        elif gt in ['1/1', '1|1']:
            AA = 0
            AB = 0
            BB = 1
        else:
            AA = 0
            AB = 1
            BB = 0

        data.append([pos, AA, AB, BB])

    data = pd.DataFrame(data)
    data.columns = ['position', 'AA', 'AB', 'BB']

    return data


def run_shapeit_grch37(
        thousand_genomes_dir,
        chromosome,
        input_bcf_file,
        output_csv,
        tempdir,
        shapeit_num_samples=100,
        shapeit_confidence_threshold=0.95
):
    meta = load_tar_meta(thousand_genomes_dir)

    snp_genotype_df = read_bcf(input_bcf_file, chromosome)

    chromosome = meta['phased_chromosome_x'] if chromosome.endswith('X') else chromosome

    if len(snp_genotype_df) == 0:
        write_null()
        return

    # Remove ambiguous positions
    snp_genotype_df = snp_genotype_df[
        (snp_genotype_df['AA'] == 1) | (snp_genotype_df['AB'] == 1) | (snp_genotype_df['BB'] == 1)]

    legend_filename = meta['legend_filename_template'].replace('{chromosome}', chromosome)
    legend_filename = os.path.join(thousand_genomes_dir, legend_filename)

    # Read snp positions from legend
    snps_df = pd.read_csv(legend_filename, compression='gzip', sep=' ',
                          usecols=['position', 'a0', 'a1'])

    # Remove indels
    snps_df = snps_df[(snps_df['a0'].isin(['A', 'C', 'T', 'G'])) & (snps_df['a1'].isin(['A', 'C', 'T', 'G']))]

    # Merge data specific inferred genotype
    snps_df = snps_df.merge(snp_genotype_df[['position', 'AA', 'AB', 'BB']], on='position', how='inner', sort=False)

    # Create genotype file required by shapeit
    snps_df['chr'] = chromosome
    snps_df['chr_pos'] = snps_df['chr'].astype(str) + ':' + snps_df['position'].astype(str)

    temp_gen_filename = os.path.join(tempdir, 'snps.gen')
    snps_df.to_csv(temp_gen_filename, sep=' ', columns=['chr', 'chr_pos', 'position', 'a0', 'a1', 'AA', 'AB', 'BB'],
                   index=False, header=False)

    temp_sample_filename = os.path.join(tempdir, 'snps.sample')
    with open(temp_sample_filename, 'w') as temp_sample_file:
        temp_sample_file.write('ID_1 ID_2 missing sex\n0 0 0 0\nUNR1 UNR1 0 2\n')

    hgraph_filename = os.path.join(tempdir, 'phased.hgraph')
    hgraph_logs_prefix = hgraph_filename + '.log'
    chr_x_flag = '--chrX' if chromosome == 'X' else ''

    genetic_map_filename = meta['genetic_map_filename_template'].replace('{chromosome}', chromosome)
    genetic_map_filename = os.path.join(thousand_genomes_dir, genetic_map_filename)

    haplotypes_filename = meta['haplotypes_filename_template'].replace('{chromosome}', chromosome)
    haplotypes_filename = os.path.join(thousand_genomes_dir, haplotypes_filename)

    sample_filename = os.path.join(thousand_genomes_dir, meta['sample_filename'])

    run_cmd(
        ['shapeit', '-M', genetic_map_filename,
         '-R', haplotypes_filename,
         legend_filename,
         sample_filename,
         '-G', temp_gen_filename,
         temp_sample_filename, '--output-graph', hgraph_filename, chr_x_flag,
         '--no-mcmc', '-L', hgraph_logs_prefix, '--seed', '12345'
         ])

    # Run shapeit to sample from phased haplotype graph
    sample_template = os.path.join(tempdir, 'sampled.{0}')
    averaged_changepoints = None

    for s in range(shapeit_num_samples):
        sample_prefix = sample_template.format(s)
        sample_log_filename = sample_prefix + '.log'
        sample_haps_filename = sample_prefix + '.haps'
        sample_sample_filename = sample_prefix + '.sample'

        success = False
        for _ in range(3):
            try:
                run_cmd(
                    ['shapeit', '-convert', '--input-graph', hgraph_filename, '--output-sample',
                     sample_prefix, '--seed', str(s), '-L', sample_log_filename]
                )
                success = True
                break
            except:
                print(f'failed sampling with seed {s}, retrying')
                continue
        if not success:
            raise Exception(f'failed to sample three times with seed {s}')

        sample_haps = pd.read_csv(sample_haps_filename, sep=' ', header=None,
                                  names=['id', 'id2', 'position', 'ref', 'alt', 'allele1', 'allele2'],
                                  usecols=['position', 'allele1', 'allele2'])
        sample_haps = sample_haps[sample_haps['allele1'] != sample_haps['allele2']]
        sample_haps['allele'] = sample_haps['allele1']
        sample_haps = sample_haps.drop(['allele1', 'allele2'], axis=1)
        sample_haps.set_index('position', inplace=True)
        sample_changepoints = sample_haps['allele'].diff().abs().astype(float).fillna(0.0)
        if averaged_changepoints is None:
            averaged_changepoints = sample_changepoints
        else:
            averaged_changepoints += sample_changepoints
        os.remove(sample_log_filename)
        os.remove(sample_haps_filename)
        os.remove(sample_sample_filename)
    averaged_changepoints /= float(shapeit_num_samples)
    last_sample_haps = sample_haps

    # Identify changepoints recurrent across samples
    changepoint_confidence = np.maximum(averaged_changepoints, 1.0 - averaged_changepoints)

    # Create a list of labels for haplotypes between recurrent changepoints
    current_hap_label = 0
    hap_label = list()

    for x in changepoint_confidence:
        if x < float(shapeit_confidence_threshold):
            current_hap_label += 1
        hap_label.append(current_hap_label)

    # Create the list of haplotypes
    haps = last_sample_haps
    haps['changepoint_confidence'] = changepoint_confidence
    haps['hap_label'] = hap_label

    haps.reset_index(inplace=True)

    haps['allele_id'] = 0

    haps_allele2 = haps.copy()
    haps_allele2['allele_id'] = 1
    haps_allele2['allele'] = 1 - haps_allele2['allele']

    haps = pd.concat([haps, haps_allele2], ignore_index=True)
    haps.sort_values(['position', 'allele_id'], inplace=True)

    haps['chromosome'] = chromosome

    haps = haps[['chromosome', 'position', 'allele', 'hap_label', 'allele_id']]

    dtypes = {
        'chromosome': str,
        'position': int,
        'allele': int,
        'hap_label': int,
        'allele_id': int
    }

    csverve.write_dataframe_to_csv_and_yaml(haps, output_csv, dtypes)


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
        run_shapeit_grch37(
            thousand_genomes_dir,
            chromosome,
            input_bcf_file,
            output_csv,
            tempdir,
            shapeit_num_samples=shapeit_num_samples,
            shapeit_confidence_threshold=shapeit_confidence_threshold
        )
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
