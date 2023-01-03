include { D4_SUM_COVERAGE                     } from '../../../modules/local/d4_sum_coverage/main'
include { D4_FILTER_COVERAGE                     } from '../../../modules/local/d4_filter_coverage/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_D4_SUM } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow CREATE_COVERAGE_MASKS {

    take:
    d4 // channel: [mandatory] [ meta, d4 ]

    main:
    ch_versions = Channel.empty()

    D4_FILTER_COVERAGE(d4)
    D4_FILTER_COVERAGE.out.bed.view()
    // TABIX_BGZIPTABIX_D4_SUM(D4_SUM_COVERAGE.out.bed)




    emit:
    bed          = D4_FILTER_COVERAGE.out.bed
    // bed      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
