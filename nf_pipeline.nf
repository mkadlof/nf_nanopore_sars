// PIPELINE PARAMETERS
// ====================
params.threads = 5                        // Number of threads (used by various modules - potentially in parallel)
params.min_cov = 50                       // Minimum coverage at a given position required to identify the variant
params.mask = 50                          // Minimum coverage at which a given position will NOT be masked in the
                                          // genome, must be less than or equal to the value specified in the -c
                                          // argument
params.bed_offset = 10                    // Value by how much we extend the amplicons in the 5' and 3' direction
params.max_depth=1200                     // The maximum coverage value considered for a position
params.length = 200                       // Minimum read length
params.mapq = 30                          // Minimal read mapping quality
params.extra_bed_offset = 10              // Value by how much we additionally extend a given amplicon in the 5' and 3'
                                          // direction,
                                          // if the required coverage was not achieved
params.quality_coverage = 1               // The minimum quality of a nucleotide required for it to be included in the
                                          // coverage count
params.pval = 0.05                        // pvalue for presence of given SNP/INDEL

// MEDAKA PARAMETERS
// ====================
params.medaka_model = 'r941_min_hac_g507' // Model used by medaka to identify SNPs, INDELs and SVs"
params.chunk_size = 800                   // Chunk size for medaka
params.chunk_overlap = 400                // Chunk overlap for medaka

params.input_genome = '/home/data/genome/sarscov2.fasta' // Path to the reference genome fasta


include { genome_id } from './modules/genome_id.nf'
include { minimap } from './modules/minimap.nf'
include { filter } from './modules/filter.nf'
include { filter as filter_2nd } from './modules/filter.nf'
include { masking } from './modules/masking.nf'
include { medaka } from './modules/medaka.nf'
include { extract_snp } from './modules/extract_snp.nf'
include { new_ref_genome } from './modules/new_ref_genome.nf'
include { coinfections } from './modules/coinfections.nf'
include { minimap_2nd } from './modules/minimap_2nd.nf'

workflow {
    // Channel
    reads = Channel.fromPath(params.reads)
    primers = Channel.value(params.primers as Path)

    // Processes
    genome_id()

    // 1st run
    minimap(reads)
    filter(minimap.out, primers)
    masking(filter.out)
    medaka(filter.out)
    extract_snp(medaka.out, genome_id.out)
    new_ref_genome(extract_snp.out)

    // 2nd run
    minimap_2nd(reads, new_ref_genome.out)
    filter_2nd(minimap_2nd.out, primers)

    // Auxiliary tasks
    coinfections(minimap.out, primers)
}
