// Declare syntax version
nextflow.enable.dsl=2



process BcftoolsConvert {
    time '24h'
    memory '2 GB'

  input:
    tuple (val(chromosome), path(vcf), path(vcf_idx))
  output:
    tuple( path("calls.vcf.gz"), path("calls.vcf.gz.csi"), val(chromosome))
 script:
    """
      bcftools view -Oz  $vcf -o calls.vcf.gz -r $chromosome
      bcftools index calls.vcf.gz
    """
}


process ShapeIt4 {
    time '24h'
    memory '2 GB'

  input:
    tuple(path(chromosome_calls_bcf), path(chromosome_calls_bcf_idx), val(chromosome), path(haplotype_reference_tar))
  output:
    val(chromosome), emit: chromosome
    path("haplotypes.csv.gz"), emit: haps_csv
    path("haplotypes.csv.gz.yaml"), emit: haps_csv_yaml
  script:
    """
        genotype_and_phase_snps infer-haps \
         --input_bcf_file $chromosome_calls_bcf \
         --thousand_genomes_tar $haplotype_reference_tar \
         --chromosome $chromosome \
         --output_csv haplotypes.csv.gz \
         --tempdir temp
    """
}


process CsverveConcat {
    time '24h'
    memory '2 GB'

  input:
    path(hap_csvs, stageAs: "?/*")
    path(hap_csv_yamls, stageAs: "?/*")
  output:
    path "merged_calls.csv.gz", emit: csv_file
    path "merged_calls.csv.gz.yaml", emit: csv_yaml_file
  script:
    """
      csverve concat --out_f merged_calls.csv.gz `echo "$hap_csvs" | sed 's/[^ ]* */--in_f &/g'`
    """
}



workflow{

    chromosomes_ch = Channel.from(params.chromosomes)
    vcf = Channel.fromPath(params.vcf)
    vcf_idx = Channel.fromPath(params.vcf_idx)
    reference_haplotyping_tar = Channel.fromPath(params.ref_tar)


    chromosomes_ch | flatten | combine(vcf)|combine(vcf_idx) | BcftoolsConvert | combine(reference_haplotyping_tar) | ShapeIt4

    CsverveConcat(ShapeIt4.out.haps_csv.collect(), ShapeIt4.out.haps_csv_yaml.collect())

    CsverveConcat.out.csv_file.subscribe { it.copyTo(params.out_csv) }
    CsverveConcat.out.csv_file.subscribe { it.copyTo(params.out_csv_yaml) }
}

