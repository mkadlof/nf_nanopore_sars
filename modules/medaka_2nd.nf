process medaka_2nd {
    publishDir "results/${sampleId}", mode: 'symlink'
    container 'ontresearch/medaka:latest'
    containerOptions "--volume ${params.genome_dir_abs_path}:/home/data/genome"

    input:
    tuple val(sampleId), path('forvariants.bam'), path('forvariants.bam.bai')
    tuple val(sampleId), path('reference_genome.fasta')

    output:
    tuple val(sampleId), path('medaka_annotated_filtered.vcf.gz')

    script:
    """
    medaka consensus --model ${params.medaka_model} \
                     --threads ${params.threads} \
                     --chunk_len ${params.chunk_size} \
                     --chunk_ovlp ${params.chunk_overlap} \
                     forvariants.bam \
                     forvariants.hdf

    medaka variant reference_genome.fasta forvariants.hdf medaka.vcf
    medaka tools annotate medaka.vcf reference_genome.fasta forvariants.bam medaka_annotated.vcf
    
    bgzip medaka.vcf; tabix medaka.vcf.gz
    bgzip medaka_annotated.vcf; tabix medaka_annotated.vcf.gz
    
    qual=`echo ${params.pval} | awk '{print 10*-log(\$1)/log(10)}'`
    bcftools filter -O z -o medaka_annotated_filtered.vcf.gz -i "GQ > \${qual} && DP >= ${params.min_cov}" medaka_annotated.vcf.gz
    """
}