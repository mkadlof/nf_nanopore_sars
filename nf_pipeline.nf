params.threads = 5

include { minimap } from './modules/minimap.nf'

workflow {
    // Channel
    reads = Channel.fromPath(params.reads)

    // Processes
    minimap(reads)
}
