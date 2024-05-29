process minimap_2nd{
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    path(reads)
    tuple val(sampleId), path('new_ref_genome_from_medaka.fasta')

    output:
    tuple val(sampleId), path('sorted_2nd.bam'), path('sorted_2nd.bam.bai')

    script:
    sampleId = reads.baseName.replace('.fastq', '')
    """
    minimap2 -a -x map-ont -t ${params.threads} new_ref_genome_from_medaka.fasta ${reads} | \
            samtools view -bS -F 2052 - | \
            samtools sort -o sorted_2nd.bam -
    samtools index sorted_2nd.bam

    """
}