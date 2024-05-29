process ambiguities {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path('forvariants.bam'), path('forvariants.bam.bai')

    output:
    tuple val(sampleId), path('detected_variants_varscan.txt')

    script:
    """
    samtools mpileup -f ${params.input_genome} -B -Q 1 forvariants.bam >> ${sampleId}.mpileup
    java -jar /opt/varscan/varscan.jar pileup2cns ${sampleId}.mpileup \
                                                  --min-avg-qual 1 \
                                                  --p-value ${params.pval} \
                                                  --min-var-freq ${params.lower_ambig} \
                                                  --min-coverage ${params.min_cov} \
                                                  --variants \
                                                  --min-reads2 0 >> detected_variants_varscan.txt
    """
}