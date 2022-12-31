//
// MOSDEPTH_CREATE_MASK
//
// Create coverage mask file for a sample set
//

include { MOSDEPTH                    } from '../../../modules/nf-core/mosdepth/main'
include { D4TOOLS_VIEW                     } from '../../../modules/local/d4tools/view/main'

workflow MOSDEPTH_CREATE_MASK {
    take:
    bam                      // channel: [mandatory] meta, bam, bai
    fasta                    // channel: [optional] fasta
    fai                      // channel: [optional] fai
    intervals                // channel: [optional]
    sample_sets              // channel: [optional]

    main:
    ch_versions = Channel.empty()

    d4 = MOSDEPTH(bam, fasta, fai).per_base_d4

    d4_intervals = d4.combine(intervals)
     	.map{ meta, bam, intervals, num_intervals ->

            // If no interval file provided (0) then add empty list
            intervals_new = num_intervals == 0 ? [] : intervals

            [[
                data_type:      meta.data_type,
                id:             meta.id,
                num_intervals:  num_intervals,
                sample:         meta.id,
            ],
             d4, intervals_new]
        }

    // This is done per interval: ideally pipe results into d4tools
    // count script that collects coverage per region -> concat
    // regions to one d4 file containing per-base coverage for sample
    // sets -> can later be used to generate mask file depending on
    // cutoff. NB: this can quickly become combinatorially expensive
    // if we want to compute coverages for every sample set; still we
    // probably want sequence masks for all sample sets (worst-case
    // scenario: even for individuals). Alternatively provide option
    // to only use population-level or study-level coverage to compute
    // cutoffs

    //D4TOOLS_MERGE(d4_intervals)

    //D4TOOLS_VIEW(d4_intervals)

    emit:
    //per_base_bed      = D4TOOLS_VIEW.out.bed           // channel: [ val(meta), [ bed ] ]
    //per_base_d4      = MOSDEPTH.out.per_base_d4           // channel: [ val(meta), [ d4 ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
