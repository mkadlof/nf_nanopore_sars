process consensus {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
     tuple val(sampleId), path('medaka_and_varscan_final_tworun.vcf.gz'), path('medaka_and_varscan_final_tworun.vcf.gz.tbi')

    output:
    tuple val(sampleId), path('medaka.fasta')

    script:
    """
    cat ${params.input_genome} | bcftools consensus medaka_and_varscan_final_tworun.vcf.gz > medaka.fasta
    """
}