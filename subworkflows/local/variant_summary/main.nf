include { VCFTOOLS as VCFTOOLS_SITE_PI } from '../../../modules/local/vcftools/main'
include { VCFTOOLS as VCFTOOLS_FREQ } from '../../../modules/local/vcftools/main'
include { VCFTOOLS as VCFTOOLS_COUNTS } from '../../../modules/local/vcftools/main'
include { VCFTOOLS as VCFTOOLS_TSTV_SUMMARY } from '../../../modules/local/vcftools/main'

include { WINDOWED_VCFTOOLS as WINDOWED_VCFTOOLS_TSTV } from '../../../modules/local/windowedvcftools/main'
include { WINDOWED_VCFTOOLS as WINDOWED_VCFTOOLS_WINDOW_PI } from '../../../modules/local/windowedvcftools/main'

workflow VARIANT_SUMMARY {

    // FIXME: better?: let mask consist of [val(meta), mask,
    // window_size], and if window_size is not "null" then run
    // windowed analyses (also?)
    take:
    vcf                // channel. [val(meta), vcf]
    mask               // channel, [val(meta), mask_bed]
    windowed_mask      // channel,  [val(meta), windowed_mask_bed]
    windows            // channel, [val(window_size)]

    main:

    windows.view()
    ch_versions = Channel.empty()

    VCFTOOLS_SITE_PI (
	vcf,
	mask
    )

    VCFTOOLS_FREQ (
	vcf,
	mask
    )

    VCFTOOLS_COUNTS (
	vcf,
	mask
    )


    VCFTOOLS_TSTV_SUMMARY (
	vcf,
	mask
    )


    // Windowed analyses
    WINDOWED_VCFTOOLS_TSTV (
	vcf,
	windowed_mask,
	windows
    )

    WINDOWED_VCFTOOLS_WINDOW_PI (
	vcf,
	windowed_mask,
	windows
    )

    // Only needed once
    ch_versions = ch_versions.mix(VCFTOOLS_SITE_PI.out.versions)

    emit:
    sites_pi      = VCFTOOLS_SITE_PI.out.sites_pi
    frq          = VCFTOOLS_FREQ.out.frq
    frq_count        = VCFTOOLS_COUNTS.out.frq_count
    tstv_summary         = VCFTOOLS_TSTV_SUMMARY.out.tstv_summary
    tstv         = WINDOWED_VCFTOOLS_TSTV.out.tstv
    windowed_pi         = WINDOWED_VCFTOOLS_WINDOW_PI.out.windowed_pi


    versions = ch_versions                     // channel: [ versions.yml ]
}
