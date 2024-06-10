process snpEff {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(consensus_vcf_gz), path(consensus_vcf_gz_tbi)

    output:
    tuple val(sampleId), path('detected_variants_consensus_annotated.vcf.gz'), path('detected_variants_consensus_annotated.txt')

    script:
    """
    java -jar /opt/snpEff/snpEff.jar ann -noStats \${GENOME_ID} \
         ${consensus_vcf_gz} > detected_variants_consensus_annotated.vcf

    bgzip --force detected_variants_consensus_annotated.vcf
    tabix detected_variants_consensus_annotated.vcf.gz

    # wywalamy pola związane z pokryciem i użyciem allelu alternation wiec zmieniamy tego awk-a
    bcftools query --format '%REF%POS%ALT| %ANN \n' \
               detected_variants_consensus_annotated.vcf.gz | \
                  cut -d "|" -f1,3,5,12 | \
                  tr "|" "\t" | \
                  awk 'BEGIN {OFS = "\t"} {if ( \$2 == "upstream_gene_variant" || \$2 == "downstream_gene_variant") {gene="."; aa="."} else {gene=\$3; aa=\$4}; print gene, \$1, aa, \$2}' > detected_variants_consensus_annotated.txt
    """
}