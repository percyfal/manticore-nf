//
// MOSDEPTH_CREATE_MASK
//

include { MOSDEPTH                    } from '../../../modules/nf-core/mosdepth/main'
include { D4TOOLS_VIEW                     } from '../../../modules/local/d4tools/view/main'

workflow MOSDEPTH_CREATE_MASK {
    take:
    bams                     // channel: [mandatory] meta, bam, bai
    fasta                    // channel: [optional] fasta
    fai                    // channel: [optional] fai
    intervals_bed_combined  // FIXME: causes MOSDEPTH fail

    main:
    ch_versions = Channel.empty()

    MOSDEPTH(bams, Channel.value([]), fasta)
    D4TOOLS_VIEW(MOSDEPTH.out.per_base_d4)

    emit:
    per_base_bed      = D4TOOLS_VIEW.out.bed           // channel: [ val(meta), [ bed ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
