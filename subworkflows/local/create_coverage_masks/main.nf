include { CAT_ZCAT as CAT_BED } from '../../../modules/local/cat/zcat/main'
include { FILTER_COVERAGE as FILTER_COVERAGE_BED } from '../../../modules/local/filter_coverage/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_COVERAGE_MASK } from '../../../modules/local/bedtools/merge/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_COVERAGE_MASK } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_PAIRS } from '../../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_SORT } from '../../../modules/nf-core/bedtools/sort/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_COVERAGE_MASK_PAIRS } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow CREATE_COVERAGE_MASKS {

    take:
    bed           // channel: [mandatory] [ meta, bed ]
    genome        // channel: path


    main:
    ch_versions = Channel.empty()

    // Filter bed files using min and max
    CAT_BED(bed)
    FILTER_COVERAGE_BED(CAT_BED.out.txt)
    BEDTOOLS_MERGE_COVERAGE_MASK(FILTER_COVERAGE_BED.out.bed)
    TABIX_BGZIPTABIX_COVERAGE_MASK(BEDTOOLS_MERGE_COVERAGE_MASK.out.bed)

    // Combine all pairs for paired analyses
    //
    // FIXME: potentially only combine coverage masks of same class
    // (auto or manual)
    //

    singlebed = TABIX_BGZIPTABIX_COVERAGE_MASK.out.gz_tbi.filter{meta, bed, tbi -> meta.sampleset.subset == true}

    pairs = singlebed
        .combine(
            singlebed
        )
        .filter{
            meta, bed, tbi, meta2, bed2, tbi2 ->
            meta.sampleset.id != meta2.sampleset.id
        }
        .map {
            meta, bed, tbi, meta2, bed2, tbi2 ->
            if (meta2.sampleset.id < meta.sampleset.id) {
                [meta2, bed2, tbi2, meta, bed, tbi]
            } else {
                [meta, bed, tbi, meta2, bed2, tbi2]
            }
        }
        .unique()

    BEDTOOLS_INTERSECT_PAIRS(
        pairs.map{
            meta, bed, tbi, meta2, bed2, tbi2 ->
            new_meta = [
                id: "${meta.id}.${meta2.id}.pair",
                mode: meta.mode,
                coverageset1: meta.coverageset,
                coverageset2: meta2.coverageset,
                sampleset1: meta.sampleset,
                sampleset2: meta2.sampleset,
            ]
            [new_meta, bed, bed2]
        }, "bed")


    BEDTOOLS_SORT(BEDTOOLS_INTERSECT_PAIRS.out.intersect, genome)
    TABIX_BGZIPTABIX_COVERAGE_MASK_PAIRS(BEDTOOLS_SORT.out.sorted)


    emit:
    gz_tbi           = TABIX_BGZIPTABIX_COVERAGE_MASK.out.gz_tbi
    pairs_gz_tbi           = TABIX_BGZIPTABIX_COVERAGE_MASK_PAIRS.out.gz_tbi

    versions = ch_versions                     // channel: [ versions.yml ]
}
