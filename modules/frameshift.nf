process frameshift {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path('medaka_and_varscan.vcf')
    val GENOME_ID

    output:
    tuple val(sampleId), path('medaka_2nd_and_varscan_final.vcf.gz'), path('medaka_2nd_and_varscan_final.vcf.gz.tbi')

    script:
    """
    min_cov_frameshift=`echo "${params.mask}*2" | bc -l`

    awk '{split(\$8,a, ";"); split(a[1], DP, "="); if(\$1 != "${GENOME_ID}" || (length(\$5) - length(\$4)) % 3 == 0) print \$0; else if (\$1 == "${GENOME_ID}" && (length(\$5) - length(\$4)) % 3 != 0 && DP[2] >= \$min_cov_frameshift) print \$0}' medaka_and_varscan.vcf >> tmp.vcf
    bcftools sort tmp.vcf | bcftools norm -c w -d all -f ${params.input_genome} | bcftools norm -c w -m -indels -f ${params.input_genome} | bcftools filter -O z -o medaka_2nd_and_varscan_final.vcf.gz -i "QUAL >= 0 && INFO/DP >= 1"
    tabix medaka_2nd_and_varscan_final.vcf.gz
    """
}