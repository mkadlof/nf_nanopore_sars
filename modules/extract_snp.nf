process extract_snp {
    publishDir "results/${sampleId}", mode: 'symlink'
    container 'ontresearch/medaka:latest'
    containerOptions "--volume ${params.genome_dir_abs_path}:/home/data/genome"

    input:
    tuple val(sampleId), path('medaka_annotated_filtered.vcf.gz')
    val REF_GENOME_ID

    output:
    tuple val(sampleId), path('medaka_annotated_filtered-2.vcf.gz')

    script:
    """
    echo "REF_GENOME_ID: ${REF_GENOME_ID}"
    zcat medaka_annotated_filtered.vcf.gz | \
        awk -v REF="${REF_GENOME_ID}" '{if(\$1 != REF || (\$1  == REF && length(\$5) == length(\$4))) print \$0}' | \
        bcftools sort | \
        bcftools norm -c w -d all -f ${params.input_genome} | \
        bcftools norm -c w -m -indels -f ${params.input_genome} | \
        bcftools filter -O z -o  medaka_annotated_filtered-2.vcf.gz -i "QUAL >= 0 && INFO/DP >= 1"

    tabix medaka_annotated_filtered-2.vcf.gz
    """
}