process extract_snp {
    publishDir "results/${sampleId}", mode: 'symlink'
//     container 'ontresearch/medaka:latest'
    containerOptions "--volume ${params.genome_dir_abs_path}:/home/data/genome"

    input:
    tuple val(sampleId), path('medaka_annotated_filtered.vcf.gz')

    output:
    tuple val(sampleId), path('medaka_and_varscan_final.vcf.gz'), path('medaka_and_varscan_final.vcf.gz.tbi')

    script:
    """
    echo "REF_GENOME_ID: >\${GENOME_ID}<"
    zcat medaka_annotated_filtered.vcf.gz | \
        awk '{if(\$1 != "\${GENOME_ID}" || (\$1 == "\${GENOME_ID}" && length(\$5) == length(\$4))) print \$0}' | \
        bcftools sort | \
        bcftools norm -c w -d all -f \${GENOME_FASTA} | \
        bcftools norm -c w -m -indels -f \${GENOME_FASTA} | \
        bcftools filter -O z -o medaka_and_varscan_final.vcf.gz -i "QUAL >= 0 && INFO/DP >= 1"

    tabix medaka_and_varscan_final.vcf.gz
    """
}