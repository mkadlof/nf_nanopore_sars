process genome_id {

    output:
    env(GENOME_ID)

    script:
    """
       GENOME_ID=`grep '>' ${params.input_genome} | head -n1 | tr -d '>' | cut -f 1 -d ' '`
    """
}