// Declare syntax version
nextflow.enable.dsl=2


process BcftoolsMpileup {
    time '24h'
    memory '2 GB'

  input:
    tuple(val(chromosome), path(input_reference), path(reference_vcf), path(input_bam), path(input_bai))
  output:
    val(chromosome), emit: chromosome
    path "chromosome_mpileup.vcf.gz", emit: vcf_file
  script:
    """
      bcftools mpileup -Oz -f $input_reference --regions-file $reference_vcf/CCDG_14151_B01_GRM_WGS_2020-08-05_${chromosome}.filtered.*.vcf.gz $input_bam -o chromosome_mpileup.vcf.gz
    """
}


process BcftoolsCall {
    time '24h'
    memory '2 GB'

  input:
    path(chromosome_mpileup_vcf)
    val(chromosome)
  output:
    tuple(val(chromosome), path("chromosome_calls.vcf.gz"), path("chromosome_calls.vcf.gz.csi"))
  script:
    """
      bcftools call -Oz -c $chromosome_mpileup_vcf -o chromosome_calls.vcf.gz
      bcftools index chromosome_calls.vcf.gz
    """
}



process ShapeIt4 {
    time '24h'
    memory '2 GB'

  input:
    tuple(val(chromosome), path(chromosome_calls_bcf), path(chromosome_calls_bcf_idx), path(haplotype_reference_dir))
  output:
    path "chromosome_phased.bcf", emit: phased_bcf
  script:
    """
        shapeit4 --input $chromosome_calls_bcf --map $haplotype_reference_dir/${chromosome}.b38.gmap.gz --region $chromosome --reference ${haplotype_reference_dir}/CCDG_14151_B01_GRM_WGS_2020-08-05_${chromosome}.filtered.*.bcf --output chromosome_phased.bcf --seed 2
    """
}


process BcftoolsConcat {
    time '24h'
    memory '2 GB'

  input:
    path(phased_bcfs, stageAs: "?/*")
  output:
    path "merged_calls.bcf", emit: bcf_file
  script:
    """
      bcftools concat -o merged_calls.bcf  $phased_bcfs
    """
}



workflow{

    chromosomes_ch = Channel.from(params.chromosomes)
    reference_ch = Channel.fromPath(params.reference)
    bam = Channel.fromPath(params.bam)
    bai = Channel.fromPath(params.bai)
    reference_haplotyping_dir = Channel.fromPath(params.ref_dir)

    chromosomes_ch | flatten | combine(reference_ch)|combine(reference_haplotyping_dir)|combine(bam)|combine(bai)| BcftoolsMpileup

    BcftoolsCall(BcftoolsMpileup.out.vcf_file, BcftoolsMpileup.out.chromosome)

    BcftoolsCall.out| combine(reference_haplotyping_dir) | ShapeIt4 | collect | BcftoolsConcat

    BcftoolsConcat.out.bcf_file.subscribe { it.copyTo(params.out_dir) }
    
}

