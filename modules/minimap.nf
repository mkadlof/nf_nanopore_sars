process minimap {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    path(reads)

    output:
    tuple val(sampleId), path('sorted.bam'), path('sorted.bam.bai')

    script:
    sampleId = reads.baseName.replace('.fastq', '')
    """
    minimap2 -a -x map-ont -t ${params.threads} -o tmp.sam /home/data/genome/sarscov2.fasta ${reads}
    samtools view -bS -F 2052 tmp.sam | samtools sort -o sorted.bam -
    samtools index sorted.bam
    """
}