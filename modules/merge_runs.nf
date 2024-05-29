process merge_runs {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path('medaka_and_varscan_final.vcf.gz'),
                         path('medaka_and_varscan_final.vcf.gz.tbi'),
                         path('medaka_2nd_and_varscan_final.vcf.gz'),
                         path('medaka_2nd_and_varscan_final.vcf.gz.tbi')

    output:
    tuple val(sampleId), path('medaka_and_varscan_final_tworun.vcf.gz'), path('medaka_and_varscan_final_tworun.vcf.gz.tbi')

    script:
    """
    bcftools concat medaka_and_varscan_final.vcf.gz medaka_2nd_and_varscan_final.vcf.gz | bcftools sort -O z -o medaka_and_varscan_final_tworun.vcf.gz
    tabix medaka_and_varscan_final_tworun.vcf.gz
    """
}