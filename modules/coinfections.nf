process coinfections {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(bam), path(bai)
    path(primers)

    output:
    tuple val(sampleId), path("${sampleId}_alternative_alleles_frequencies.png") , path("${sampleId}_coinfections_summary.txt")

    script:
    """
    length_filtering=`echo "${params.length} - 40 " | bc -l`

    # We take the original unfiltered BAM, but we mask the primers.
    samtools ampliconclip --filter-len \${length_filtering} --both-ends -b ${primers} --tolerance ${params.bed_offset} -o forcoinfections_presorted.bam -O bam sorted.bam
    samtools sort  -@ ${params.threads} -o forcoinfections.bam forcoinfections_presorted.bam
    samtools index forcoinfections.bam
    samtools mpileup  -f \${GENOME_FASTA} -Q 1 forcoinfections.bam >> coinfections.mpileup

    # Calling varscan
    java -jar /opt/varscan/varscan.jar pileup2snp coinfections.mpileup  --min-avg-qual 1 --p-value 0.5 --min-var-freq 0.05 --min-coverage ${params.min_cov} --min-reads2 0  >> detected_variants_varscan_coinfections.txt

    # Script creates results based on txt file and template file from EQA2023 results
    predict_coinfections_nanopore.py detected_variants_varscan_coinfections.txt ${sampleId} /home/data/coinfections/ESIB_EQA_2023.SARS1.04_coinfections.txt /home/data/coinfections/ESIB_EQA_2023.SARS1.05_coinfections.txt
    """
}