//
// MOSDEPTH_COVERAGE
//
// Summarize coverage with mosdepth and d4tools for all samplesets.
// NB! mosdepth per-base coverage *always* outputs the entire genome,
// so we cannot partition on intervals here.
//

include { MOSDEPTH                    } from '../../../modules/nf-core/mosdepth/main'
include { BEDTOOLS_UNIONBEDG          } from '../../../modules/local/bedtools/unionbedg/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_UNIONBEDG } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow MOSDEPTH_COVERAGE {
    take:
    bam                      // channel: [mandatory] meta, bam, bai
    samplesets               // channel: [mandatory] meta, samples
    fasta                    // channel: [optional] fasta
    fai                      // channel: [optional] fai
    intervals                // channel: [optional] intervals


    main:
    ch_versions = Channel.empty()


    MOSDEPTH(bam, intervals, fasta)
    // Convert samplesets to length 12: map mosdepth output to join to
    // samplesets output and then group by tuple
    bed = MOSDEPTH.out.per_base_bed.map{
	[it[0].id, it[1]]
    }.cross(
	samplesets.map{it -> it[1].collect{key -> [key, it[0].id]}}.flatten().collate(2)
    ).map{[[id:it[1][1]], it[1][0], it[0][1]]}
	.groupTuple()
	.map{
	    [[id:it[0].id, samples:it[1]], it[2]]
	}
    // Possibly do this by region
    BEDTOOLS_UNIONBEDG(bed, fai)
    TABIX_BGZIPTABIX_UNIONBEDG(BEDTOOLS_UNIONBEDG.out.bed)

    emit:
    gz_tbi       = TABIX_BGZIPTABIX_UNIONBEDG.out.gz_tbi
    versions = ch_versions                     // channel: [ versions.yml ]
}
