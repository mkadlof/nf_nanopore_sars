process minimap {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    path(reads)

    output:
    tuple val(sampleId), path('sorted.bam'), path('sorted.bam.bai')

    script:
    sampleId = reads.baseName.replace('.fastq', '')
    """
    minimap2 -a -x map-ont -t ${params.threads} ${params.input_genome} ${reads} | \
        samtools view -bS -F 2052 - | \
        samtools sort -o sorted.bam -
    samtools index sorted.bam
    """
}