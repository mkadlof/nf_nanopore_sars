process simpleStats {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path("consensus_masked.fa")

    output:
    tuple val(sampleId), path("N_summary.txt")

    script:
    """
    calculate_N.py consensus_masked.fa N_summary.txt
    """
}