process wgsMetrics {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path("forvariants.bam"), path("forvariants.bam.bai")

    output:
    tuple val(sampleId), path("picard_statistics.txt")

    script:
    """
    java -jar /opt/picard/picard.jar CollectWgsMetrics \
                  -I forvariants.bam -R \${GENOME_FASTA} \
                  -O picard_statistics.txt \
                  -Q ${params.quality_coverage} \
                  -MQ ${params.mapq} \
                  -COUNT_UNPAIRED TRUE
    """
}