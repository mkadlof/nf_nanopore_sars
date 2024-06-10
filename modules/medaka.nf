process medaka {
    publishDir "results/${sampleId}", mode: 'symlink'
    container 'ontresearch/medaka:latest'
    containerOptions "--volume ${params.genome_dir_abs_path}:/home/data/genome"

    input:
    tuple val(sampleId), path('forvariants.bam'), path('forvariants.bam.bai')

    output:
    tuple val(sampleId), path('medaka_annotated_filtered.vcf.gz')

    script:
    """
    GENOME_ID="MN908947.3"
    GENOME_FASTA="/home/data/genome/sarscov2.fasta"
    medaka consensus --model ${params.medaka_model} \
                     --threads ${params.threads} \
                     --chunk_len ${params.chunk_size} \
                     --chunk_ovlp ${params.chunk_overlap} \
                     forvariants.bam \
                     forvariants.hdf

    medaka variant \${GENOME_FASTA} forvariants.hdf medaka.vcf
    medaka tools annotate medaka.vcf \${GENOME_FASTA} forvariants.bam medaka_annotated.vcf
    
    bgzip medaka.vcf; tabix medaka.vcf.gz
    bgzip medaka_annotated.vcf; tabix medaka_annotated.vcf.gz
    
    qual=`echo ${params.pval} | awk '{print 10*-log(\$1)/log(10)}'`
    bcftools filter -O z -o medaka_annotated_filtered.vcf.gz -i "GQ > \${qual} && DP >= ${params.min_cov}" medaka_annotated.vcf.gz
    """
}