process kraken2 {
    publishDir "results/${sampleId}", mode: 'symlink'
    containerOptions "--volume ${params.kraken2_db_absolute_path_on_host}:/home/external_databases/kraken2"
    maxForks 5

    input:
    path(reads)

    output:
    tuple val(sampleId), path('raport_kraken2.txt'), path('raport_kraken2_individualreads.txt')

    script:
    sampleId = reads.baseName.replace('.fastq', '')
    """
    kraken2 --db /home/external_databases/kraken2 \
            --report raport_kraken2.txt \
            --threads ${params.threads} \
            --gzip-compressed \
            --minimum-base-quality 30 \
            --use-names ${reads} >> raport_kraken2_individualreads.txt
    """
}