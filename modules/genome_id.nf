process genome_id {

    output:
    stdout

    script:
    """
        grep '>' ${params.input_genome} | head -n1 | tr -d '>' | cut -f 1 -d ' '
    """
}