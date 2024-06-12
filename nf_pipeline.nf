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
params.pval = 0.05                        // p-value for presence of given SNP/INDEL
params.lower_ambig = 0.45                 // A value in the range 0-1, the percentage of reads with an alternative
                                          // allele version at a given position required to introduce the ambiguous
                                          // symbol at that position.
params.upper_ambig = 0.55                 // A value in the range 0-1, the percentage of reads with an alternative
                                          // allele version at a given position required to introduce the alternative
                                          // allele at that position, and not the ambiguous allele symbol. The value
                                          // must be greater than that provided with the lower_ambig argument.


// MEDAKA PARAMETERS
// ====================
params.medaka_model = 'r941_min_hac_g507' // Model used by medaka to identify SNPs, INDELs and SVs"
params.chunk_size = 800                   // Chunk size for medaka
params.chunk_overlap = 400                // Chunk overlap for medaka


include { minimap } from './modules/minimap.nf'
include { filter } from './modules/filter.nf'
include { filter as filter_2nd } from './modules/filter.nf'
include { lowCov } from './modules/low_cov.nf'
include { medaka } from './modules/medaka.nf'
include { medaka_2nd } from './modules/medaka_2nd.nf'
include { extract_snp } from './modules/extract_snp.nf'
include { new_ref_genome } from './modules/new_ref_genome.nf'
include { coinfections } from './modules/coinfections.nf'
include { minimap_2nd } from './modules/minimap_2nd.nf'
include { ambiguities } from './modules/ambiguities.nf'
include { merge } from './modules/merge.nf'
include { frameshift } from './modules/frameshift.nf'
include { merge_runs } from './modules/merge_runs.nf'
include { consensus } from './modules/consensus.nf'
include { consensusMasking } from './modules/consensusMasking.nf'
include { pangolin } from './modules/pangolin.nf'
include { nextclade } from './modules/nextclade.nf'
include { snpEff } from './modules/snpEff.nf'
include { kraken2 } from './modules/kraken2.nf'
include { modeller } from './modules/modeller.nf'
include { dehumanization } from './modules/dehumanization.nf'
include { simpleStats } from './modules/simpleStats.nf'

workflow {
    // Channel
    reads = Channel.fromPath(params.reads)
    primers = Channel.value(params.primers as Path)

    // Processes

    // 1st run
    minimap(reads)
    filter(minimap.out, primers)
    lowCov(filter.out)
    medaka(filter.out)
    extract_snp(medaka.out)
    new_ref_genome(extract_snp.out)

    // 2nd run
    minimap_2nd(reads, new_ref_genome.out)
    filter_2nd(minimap_2nd.out, primers)
    medaka_2nd(filter_2nd.out, new_ref_genome.out)

    // post-runs tasks
    ambiguities(filter_2nd.out)
    merge(medaka_2nd.out.join(ambiguities.out))
    frameshift(merge.out)
    merge_runs(extract_snp.out.join(frameshift.out))
    consensus(merge_runs.out)
    consensusMasking(lowCov.out[1].join(consensus.out))

    // Variant identification
    pangolin(consensusMasking.out)
    nextclade(consensusMasking.out)

    // Auxiliary tasks
    kraken2(reads)
    dehumanization(minimap.out, reads)
    simpleStats(consensusMasking.out)

    coinfections(minimap.out, primers)
    snpEff(merge_runs.out.join(ambiguities.out))
    modeller(nextclade.out[1])
}
