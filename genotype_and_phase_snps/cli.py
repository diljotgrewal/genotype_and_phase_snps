"""Console script for csverve."""

import click

import genotype_and_phase_snps.utils as utils
import csverve.api as csverve

@click.group()
def cli():
    pass


@cli.command()
@click.option('--input_bcf_file', required=True, help='CSV file path, allows multiple paths.')
@click.option('--output_csv', required=True, help='Path of resulting merged CSV.')
@click.option('--thousand_genomes_tar', required=True, help='How to join CSVs.')
@click.option('--chromosome', required=True, help='Column to join CSVs on, allowes multiple.')
@click.option('--tempdir', required=True, help='Writer header to resulting CSV.')
@click.option('--shapeit_num_samples', default=100, help='Writer header to resulting CSV.')
@click.option('--shapeit_confidence_threshold', default=0.95, help='Writer header to resulting CSV.')
def infer_haps(
        input_bcf_file,
        output_csv,
        thousand_genomes_tar,
        chromosome,
        tempdir,
        shapeit_num_samples=100,
        shapeit_confidence_threshold=0.95
):
    utils.infer_haps(
        input_bcf_file,
        output_csv,
        thousand_genomes_tar,
        chromosome,
        tempdir,
        shapeit_num_samples=shapeit_num_samples,
        shapeit_confidence_threshold=shapeit_confidence_threshold
    )

if __name__ == "__main__":
    cli()
