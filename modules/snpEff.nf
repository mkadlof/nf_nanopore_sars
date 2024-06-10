process snpEff {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(consensus_vcf_gz), path(consensus_vcf_gz_tbi)

    output:
    tuple val(sampleId), path('detected_variants_consensus_annotated.vcf.gz'), path('detected_variants_medaka_annotated.txt')

    script:
    """
    echo a
    java -jar /opt/snpEff/snpEff.jar ann -noStats \${GENOME_ID} ${consensus_vcf_gz} > detected_variants_consensus_annotated.vcf
    echo b

    bgzip --force detected_variants_consensus_annotated.vcf
    tabix detected_variants_consensus_annotated.vcf.gz

    ## Niestety po ominięciu longshot-a nie mamy juz dostępu do użycia alleli z pliku vcf wiec przygotowanie pliku txt
    ## jest 3 etapowe.
    ## Part1 wyciągamy gen, mutacje, efekt białkowy i pokrycie na danej pozycji z vcf-a.
    bcftools query -f '%POS | %REF%POS%ALT | %DP | %ANN \n' detected_variants_medaka_annotated.vcf.gz |\
                    tr "|" "\t" |\
                    cut -f1,2,3,5,7,14 |\
                    awk 'BEGIN {OFS = "\t"}  {if ( \$4 == "upstream_gene_variant" || \$4 == "downstream_gene_variant") {gene="."; aa="."} else {gene=\$5; aa=\$6};  print \$1, gene, \$2, aa, \$3}' |\
                    sort -k1  >> part1.txt

    cat detected_variants_varscan.txt  | grep -v Pos |  cut -f2,7 | sort -k1 >> part2.txt

    ## Łączenie to jest przedziwna konstrukcja przy pomocy join-a
    ## -1 1 i -2 1 oznaczają ze łączymy po pierwszej kolumnie w plikach 1 i 2 -o 1.1,2.1 to output gdzie w konstrukcji
    ## 1.3 "1" oznacza numer pliku a "3" kolumnę w pliku     ## -a 1 oznacza by wypisywał również linie z pliku 1,
    ## które nie maja klucza w pliku 2, a  "-e" mówi jaki znak wtedy wypisać robię to gdy mutacja z medaki nie było
    ## w varscan

    join -a 1 -1 1 -2 1 -o1.1,1.2,1.3,1.4,1.5,2.2  -e '-'  part1.txt part2.txt | sort -nk1 | cut -d " " -f2- | tr " " "\t" >> detected_variants_medaka_annotated.txt

    rm part1.txt part2.txt
    """
}