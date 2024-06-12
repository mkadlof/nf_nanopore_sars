process dehumanization  {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai')
    path(reads)

    output:
    tuple val(sampleId), path('reads_nohuman.fq.gz')

    script:
    """
    samtools view mapped_reads.bam | cut -f1 | sort | uniq >> lista_id_nohuman.txt
    seqtk subseq ${reads} lista_id_nohuman.txt >> reads_nohuman.fq

    gzip reads_nohuman.fq
    """
}