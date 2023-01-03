/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowManticore.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.multiqc_config,
    params.vcf,
    params.fasta,
    params.fasta_fai,
    params.intervals,
    params.dict,
    params.roi_fof,
    params.sample_sets
]

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.intervals && !params.intervals.endsWith("bed") && !params.intervals.endsWith("list")) exit 1, "Intervals file must end with .bed, .list, or .interval_list"

// Process roi_fof file
if (params.roi_fof) {
    ch_roi = Channel.from(file(params.roi_fof)).splitCsv(header: false, strip: true)
	.map{ row ->
	    if (!(row[0] == "site" || row[0] == "window")) {
		log.error("first column ('${row[0]}') in roi_fof ($params.roi_fof) must be either 'site' or 'window'")
		System.exit(1)
	    }
	    if (!row[1].endsWith("bed")) exit 1, "in ${params.roi_fof}: ${row[1]}: Intervals must be in bed format"
	    if (!file(row[1], checkIfExists: true)) exit 1, "no such file ${row[1]}"
	    [[id:file(row[1]).baseName, mode: row[0]], file(row[1])]
	}
} else {
    ch_roi = Channel.empty()
}
// Process sample sets file
if (params.sample_sets) {
    ch_sample_sets = Channel.from(file(params.sample_sets)).splitCsv(sep: "\t", header: false, strip: true)
	.map{ row ->
	    if (row[0] == "ALL") exit 1, "in ${params.sample_sets}: 'ALL' is a reserved keyword!"
	    [[id: row[0]], row[1].split(",")]
	}
} else {
    ch_sample_sets = Channel.empty()
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

fasta              = params.fasta              ? Channel.fromPath(params.fasta).collect()                    : Channel.empty()
fasta_fai          = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).collect()                : Channel.empty()
vcf                = params.vcf                ? Channel.fromPath(params.vcf).collect()                      : Channel.empty()
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK                                  } from '../subworkflows/local/input_check'
include { PREPARE_GENOME                               } from '../subworkflows/local/prepare_genome/main'
include { PREPARE_INTERVALS                            } from '../subworkflows/local/prepare_intervals/main'
include { MOSDEPTH_COVERAGE                            } from '../subworkflows/local/mosdepth_coverage/main'
include { CREATE_COVERAGE_MASKS                        } from '../subworkflows/local/create_coverage_masks/main'
include { CREATE_MASKS                                 } from '../subworkflows/local/create_masks/main'
include { VARIANT_SUMMARY                              } from '../subworkflows/local/variant_summary/main'
include { D4UTILS_SUM                                  } from '../modules/local/d4utils/sum/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MANTICORE {

    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    // Build indices
    PREPARE_GENOME(
        fasta,
        fasta_fai
    )

    // Gather built indices and resources
    dict               = params.fasta                   ? params.dict                       ? Channel.fromPath(params.dict).collect()                  : PREPARE_GENOME.out.dict                  : []
    fasta_fai          = params.fasta                   ? params.fasta_fai                  ? Channel.fromPath(params.fasta_fai).collect()             : PREPARE_GENOME.out.fasta_fai             : []

    // Build intervals if needed
    PREPARE_INTERVALS(fasta_fai)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    intervals_bed_combined      = params.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_combined  // [interval.bed] all intervals in one file
    intervals                   = PREPARE_INTERVALS.out.intervals_bed        // [interval, num_intervals] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi        = PREPARE_INTERVALS.out.intervals_bed_gz_tbi // [interval_bed, tbi, num_intervals] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather

    // Gather used softwares versions
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
    ch_versions = ch_versions.mix(PREPARE_INTERVALS.out.versions)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    ch_sample_sets_all = ch_sample_sets.mix(INPUT_CHECK.out.bams.map{[[id: "ALL"], it[0].id]}.groupTuple())

    MOSDEPTH_COVERAGE(
        INPUT_CHECK.out.bams,
        fasta,
        fasta_fai,
        intervals
    )
    ch_versions = ch_versions.mix(MOSDEPTH_COVERAGE.out.versions.first())

    // Apply (multiple) coverage cutoffs to each and every sampleset;
    // need to make all combinations!
    d4 = MOSDEPTH_COVERAGE.out.d4
    D4UTILS_SUM(d4)
    CREATE_COVERAGE_MASKS(
	ch_sample_sets_all.combine(d4).map{
	    [[id: "${it[0].id}-${it[2].id}", min: null, max: null, samples: it[1], interval: it[2].id], it[3]]}.
	groupTuple()
    )
    //
    // SUBWORKFLOW: CREATE_MASKS: create mask files from input mask and
    // coverages
    //
    // FIXME: Combine masks with coverage cutoffs for all masks; will
    // generate lots of large genome mask files...
    CREATE_MASKS (
	PREPARE_GENOME.out.genome_bed, fasta,
	ch_roi,
    )
    ch_versions = ch_versions.mix(CREATE_MASKS.out.versions.first())

    // Expand mask file list to modes and combine with window sizes
    window_sizes = Channel.of(params.window_sizes.split(","))
    masks = CREATE_MASKS.out.fasta.map{
	[[id: "${it[0].id}", bed: "${it[0].bed}"],
	 it[0].mode,
	 it[1]]
    }.transpose(by: 1).map{
	[[id: it[0].id, bed: it[0].bed, mode: it[1]], it[2]]
    }.combine(window_sizes).map{
	[[id: it[0].id, bed: it[0].bed, mode: it[0].mode, window_size: it[0].mode == "window" ? it[2] : null], it[1]]
    }.unique()

    VARIANT_SUMMARY(
	vcf.map{[[id:"nucleotide_diversity.${vcf.baseName}"], it]},
	masks
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowManticore.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowManticore.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
