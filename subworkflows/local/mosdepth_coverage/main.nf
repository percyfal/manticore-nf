//
// MOSDEPTH_COVERAGE
//
// Summarize coverage with mosdepth and d4tools
//

include { MOSDEPTH                    } from '../../../modules/nf-core/mosdepth/main'
include { D4TOOLS_MERGE                     } from '../../../modules/local/d4tools/merge/main'
include { D4_SUM_COVERAGE                     } from '../../../modules/local/d4_sum_coverage/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_D4_SUM } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow MOSDEPTH_COVERAGE {
    take:
    bam                      // channel: [mandatory] meta, bam, bai
    fasta                    // channel: [optional] fasta
    fai                      // channel: [optional] fai
    intervals

    main:
    ch_versions = Channel.empty()

    bam_intervals = bam.combine(intervals)
	.map{ meta, bam, bai, intervals, num_intervals ->
	    intervals_new = num_intervals == 0 ? [] : intervals

	    [intervals_new]
	}

    mosdepth_bams = bam.combine(intervals)
	.map {meta, bam, bai, intervals, num_intervals ->
	    basename = intervals.name
	    [[
	     	data_type: meta.data_type,
	     	id: "$meta.id-$basename", // Must be unique!
	     	sample: meta.id,
		interval: basename,
		num_intervals: num_intervals,
	    ], bam, bai]
	}

    MOSDEPTH(mosdepth_bams, bam_intervals, fasta).per_base_d4

    MOSDEPTH.out.per_base_d4.branch{
	intervals:    it[0].num_intervals > 1
	no_intervals: it[0].num_intervals <= 1

    }.set{ mosdepth_d4_branch }

    // FIXME: new_meta should contain label for sample set
    D4TOOLS_MERGE(
	mosdepth_d4_branch.intervals
	    .map{ meta, d4 ->
		new_meta = [
		    id: meta.interval,
		]
		label = "$d4:$meta.sample"
		[new_meta, d4, label]
	    }.groupTuple()
    )

    D4_SUM_COVERAGE(D4TOOLS_MERGE.out.d4)

    TABIX_BGZIPTABIX_D4_SUM(D4_SUM_COVERAGE.out.bed)


    emit:
    per_base_d4      = D4TOOLS_MERGE.out.d4           // channel: [ val(meta), [ d4 ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
