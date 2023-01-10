include { VCFTOOLS as VCFTOOLS_SITE_PI } from '../../../modules/local/vcftools/main'
include { VCFTOOLS as VCFTOOLS_FREQ } from '../../../modules/local/vcftools/main'
include { VCFTOOLS as VCFTOOLS_COUNTS } from '../../../modules/local/vcftools/main'
include { VCFTOOLS as VCFTOOLS_TSTV_SUMMARY } from '../../../modules/local/vcftools/main'
include { VCFTOOLS_FST as MAIN_VCFTOOLS_FST} from '../../../modules/local/vcftools_fst/main'

include { WINDOWED_VCFTOOLS as WINDOWED_VCFTOOLS_TSTV } from '../../../modules/local/windowedvcftools/main'
include { WINDOWED_VCFTOOLS as WINDOWED_VCFTOOLS_WINDOW_PI } from '../../../modules/local/windowedvcftools/main'
include { WINDOWED_VCFTOOLS as WINDOWED_VCFTOOLS_TAJIMAS } from '../../../modules/local/windowedvcftools/main'
include { VCFTOOLS_FST as WINDOWED_VCFTOOLS_FST } from '../../../modules/local/vcftools_fst/main'


workflow VARIANT_SUMMARY {

    take:
    vcf                // channel. [val(meta), vcf]
    mask               // channel, [val(meta), mask_fasta]
    pairs              // channel, [val(meta), pairs]

    main:

    ch_versions = Channel.empty()

    mask.branch{
        site: it[0].mode == "site"
        window: it[0].mode == "window"
    }.set{mask_mode}

    VCFTOOLS_SITE_PI (
        vcf,
        mask_mode.site
    )

    VCFTOOLS_FREQ (
        vcf,
        mask_mode.site
    )

    VCFTOOLS_COUNTS (
        vcf,
        mask_mode.site
    )

    VCFTOOLS_TSTV_SUMMARY (
        vcf,
        mask_mode.site
    )

    WINDOWED_VCFTOOLS_TSTV (
        vcf,
        mask_mode.window
    )

    WINDOWED_VCFTOOLS_WINDOW_PI (
        vcf,
        mask_mode.window
    )

    WINDOWED_VCFTOOLS_TAJIMAS (
        vcf,
        mask_mode.window

    )

    // Pairs analyses
    pairs.branch{
        site: it[0].mode == "site"
        window: it[0].mode == "window"
    }.set{pairs_mode}

    MAIN_VCFTOOLS_FST(
        vcf,
        pairs_mode.site
    )

    WINDOWED_VCFTOOLS_FST(
        vcf,
        pairs_mode.window
    )

    // Only needed once
    ch_versions = ch_versions.mix(VCFTOOLS_SITE_PI.out.versions)

    emit:
    sites_pi      = VCFTOOLS_SITE_PI.out.sites_pi
    frq          = VCFTOOLS_FREQ.out.frq
    frq_count        = VCFTOOLS_COUNTS.out.frq_count
    tstv_summary         = VCFTOOLS_TSTV_SUMMARY.out.tstv_summary
    fst            = MAIN_VCFTOOLS_FST.out.weir_fst
    tstv         = WINDOWED_VCFTOOLS_TSTV.out.tstv
    windowed_pi         = WINDOWED_VCFTOOLS_WINDOW_PI.out.windowed_pi
    windowed_fst        = WINDOWED_VCFTOOLS_FST.out.weir_fst

    versions = ch_versions                     // channel: [ versions.yml ]
}
