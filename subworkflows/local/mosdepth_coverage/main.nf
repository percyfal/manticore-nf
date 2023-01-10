//
// MOSDEPTH_COVERAGE
//
// Summarize coverage with mosdepth for all samplesets. NB! mosdepth
// per-base coverage *always* outputs the entire genome, so we cannot
// partition on intervals here.
//

include { MOSDEPTH                    } from '../../../modules/nf-core/mosdepth/main'
include { BEDTOOLS_UNIONBEDG          } from '../../../modules/local/bedtools/unionbedg/main'
include { SUM_COVERAGE          } from '../../../modules/local/sum_coverage/main'
include { COVERAGE_STATS          } from '../../../modules/local/coverage_stats/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_UNIONBEDG } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_SUM_COVERAGE } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow MOSDEPTH_COVERAGE {
    take:
    bam                      // channel: [mandatory] meta, bam, bai
    samplesets               // channel: [mandatory] meta, samples
    fasta                    // channel: [optional] fasta
    fai                      // channel: [optional] fai
    genome_txt               // channel: [optional] txt
    intervals                // channel: [optional] intervals

    main:
    ch_versions = Channel.empty()


    MOSDEPTH(bam, intervals, fasta)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    sample2set = samplesets.map{meta -> meta.samples.collect{key -> [meta.id, key]}}
        .flatten()
        .collate(2)
        .combine(
            samplesets.map{meta -> [meta.id, meta]}, by:0
        )
        .map{
            id, sample, sampleset ->
            [sample, sampleset]
        }

    bed = MOSDEPTH.out.per_base_bed
        .map{
            meta, bed ->
            [meta.id, bed]
        }
        .cross(sample2set)
        .map{
            bed, sample ->
            [sample[1], bed[1]]
        }
        .groupTuple()

    // FIXME: add possibility to do this by interval scatter-gather
    BEDTOOLS_UNIONBEDG(bed, fai)
    ch_versions = ch_versions.mix(BEDTOOLS_UNIONBEDG.out.versions)
    SUM_COVERAGE(BEDTOOLS_UNIONBEDG.out.bed)
    COVERAGE_STATS(SUM_COVERAGE.out.bed)
    ch_versions = ch_versions.mix(SUM_COVERAGE.out.versions)
    ch_versions = ch_versions.mix(COVERAGE_STATS.out.versions)

    TABIX_BGZIPTABIX_UNIONBEDG(
        BEDTOOLS_UNIONBEDG.out.bed.map{
            sampleset, bed ->
            [[id:"${sampleset.id}.unionbedg", sampleset: sampleset], bed]
        })

    TABIX_BGZIPTABIX_SUM_COVERAGE(SUM_COVERAGE.out.bed)
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_UNIONBEDG.out.versions)
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX_SUM_COVERAGE.out.versions)

    gz_tbi = TABIX_BGZIPTABIX_UNIONBEDG.out.gz_tbi
    // For some reason sampleset_id gets converted to sampleset unless
    // we apply this mapping?!?

    sum_gz_tbi = TABIX_BGZIPTABIX_SUM_COVERAGE.out.gz_tbi
        .map{
            sampleset, bed, gz ->
            [[id: "${sampleset.id}", sampleset: sampleset], bed]
        }

    coverage_stats = COVERAGE_STATS.out.stats.map{
        meta, x ->
        new_meta = [
            id: "${meta.id}.auto",
            coverage_id: "auto",
            sampleset: meta,
            nsites: x.split("\t")[0].toInteger(),
            mean: x.split("\t")[1].toFloat(),
            var: x.split("\t")[2].toFloat(),
            sd: x.split("\t")[3].toFloat(),
        ]
    }

    emit:
    gz_tbi           = gz_tbi
    sum_gz_tbi       = sum_gz_tbi
    stats_txt        = COVERAGE_STATS.out.txt
    stats        = coverage_stats
    versions = ch_versions                     // channel: [ versions.yml ]
}
