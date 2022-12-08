// Declare syntax version
nextflow.enable.dsl=2



process BcftoolsConvert {
    time '24h'
    memory '2 GB'

  input:
    tuple (val(chromosome), path(vcf), path(vcf_idx))
  output:
    tuple( path("calls.bcf"), path("calls.bcf.csi"), val(chromosome))
 script:
    """
      bcftools view -Ob  $vcf -o calls.bcf -r $chromosome
      bcftools index calls.bcf
    """
}


process ShapeIt4 {
    time '24h'
    memory '2 GB'

  input:
    tuple(path(chromosome_calls_bcf), path(chromosome_calls_bcf_idx), val(chromosome), path(haplotype_reference_dir))
  output:
    val(chromosome), emit: chromosome
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
    vcf = Channel.fromPath(params.vcf)
    vcf_idx = Channel.fromPath(params.vcf_idx)
    reference_haplotyping_dir = Channel.fromPath(params.ref_dir)

    //chromosomes_ch.combine(vcf).combine(vcf_idx).view()
    //BcftoolsConvert(tuple(chromosomes_ch, vcf, vcf_idx))


    chromosomes_ch | flatten | combine(vcf)|combine(vcf_idx) | BcftoolsConvert | combine(reference_haplotyping_dir) | ShapeIt4

    BcftoolsConcat(ShapeIt4.out.phased_bcf.collect())   
   
    BcftoolsConcat.out.bcf_file.subscribe { it.copyTo(params.out_dir) }
}

