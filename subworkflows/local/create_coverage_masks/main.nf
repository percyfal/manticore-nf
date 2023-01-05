include { CAT_ZCAT as CAT_BED } from '../../../modules/local/cat/zcat/main'
include { FILTER_COVERAGE as FILTER_COVERAGE_BED } from '../../../modules/local/filter_coverage/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_COVERAGE_MASK } from '../../../modules/local/bedtools/merge/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_COVERAGE_MASK } from '../../../modules/nf-core/tabix/bgziptabix/main'


workflow CREATE_COVERAGE_MASKS {

    take:
    bed           // channel: [mandatory] [ meta, bed ]


    main:
    ch_versions = Channel.empty()

    // Filter bed files using min and max
    CAT_BED(bed)
    FILTER_COVERAGE_BED(CAT_BED.out.txt)
    BEDTOOLS_MERGE_COVERAGE_MASK(FILTER_COVERAGE_BED.out.bed)
    TABIX_BGZIPTABIX_COVERAGE_MASK(BEDTOOLS_MERGE_COVERAGE_MASK.out.bed)

    emit:
    gz_tbi           = TABIX_BGZIPTABIX_COVERAGE_MASK.out.gz_tbi

    versions = ch_versions                     // channel: [ versions.yml ]
}
