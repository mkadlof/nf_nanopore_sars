process new_ref_genome {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path('medaka_annotated_filtered-2.vcf.gz'), path('medaka_annotated_filtered-2.vcf.gz.tbi')

    output:
    tuple val(sampleId), path('new_ref_genome_from_medaka.fasta')

    script:
    """
    cat ${params.input_genome} | bcftools consensus medaka_annotated_filtered-2.vcf.gz  > new_ref_genome_from_medaka.fasta
    """
}