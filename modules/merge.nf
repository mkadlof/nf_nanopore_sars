process merge {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path('medaka_2nd_annotated_filtered.vcf.gz'), path('detected_variants_varscan.txt')

    output:
    tuple val(sampleId), path('medaka_and_varscan.vcf')

    script:
    """
    merge_varscan_with_medaka_final.py medaka_2nd_annotated_filtered.vcf.gz \
                                       detected_variants_varscan.txt \
                                       ${params.upper_ambig} \
                                       medaka_and_varscan.vcf \
                                       ${params.min_cov}
    """
}