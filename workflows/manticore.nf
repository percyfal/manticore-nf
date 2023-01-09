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
    params.sample_sets,
]

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

Channel.fromPath(params.input).splitCsv(header: true, sep: ',')
    .map{ row ->
        WorkflowManticore.validateSampleRow(workflow, row, log)
        create_bam_channel(row)
    }
    .set{ch_input}

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
    ch_sample_sets_input = Channel.from(file(params.sample_sets)).splitCsv(sep: "\t", header: false, strip: true)
        .map{ row ->
            if (row[0] == "ALL") exit 1, "in ${params.sample_sets}: 'ALL' is a reserved keyword!"
            [[id: row[0], subset: true], row[1].split(",")]
        }
} else {
    ch_sample_sets_input = Channel.empty()
}
// Add ALL sample set
ch_sample_sets = ch_sample_sets_input.mix(ch_input.map{[[id: "ALL", subset: false], it[0].id]}.groupTuple())

// Manual coverage parameter sets
if (params.coverage_min) {
    coverage_min = Channel.from(params.coverage_min.split(",")).map{[[id:it.split(":")[0]], [min: it.split(":")[1].toFloat()]]}
} else {
    coverage_min = Channel.empty()
}
if (params.coverage_max) {
    coverage_max = Channel.from(params.coverage_max.split(",")).map{[[id:it.split(":")[0]], [max: it.split(":")[1].toFloat()]]}
} else {
    coverage_max = Channel.empty()
}

ch_coverage_set = coverage_min.join(coverage_max, remainder: true)
    .join(ch_sample_sets.map{meta, samples -> [[id: meta.id], meta.subset, samples]})
    .map{
        meta, min, max, subset, samples ->
        new_meta = [
            id: "${meta.id}.manual",
            subset: subset,
            coverage_id: "manual",
            sampleset_id: meta.id,
            samples: samples,
            min: min? min.min:null,
            max: max? max.max:null,
        ]
    }
    .mix(ch_sample_sets.map{
            meta, samples ->
            new_meta = [
                id: "${meta.id}.auto",
                subset: meta.subset,
                coverage_id: "auto",
                sampleset_id: meta.id,
                samples: samples,
                min: null,
                max: null,
            ]
        }
    )

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

include { PREPARE_GENOME                               } from '../subworkflows/local/prepare_genome/main'
include { PREPARE_INTERVALS                            } from '../subworkflows/local/prepare_intervals/main'
include { MOSDEPTH_COVERAGE                            } from '../subworkflows/local/mosdepth_coverage/main'
include { CREATE_COVERAGE_MASKS                        } from '../subworkflows/local/create_coverage_masks/main'
include { CREATE_SEQUENCE_MASKS                                 } from '../subworkflows/local/create_sequence_masks/main'
include { VARIANT_SUMMARY                              } from '../subworkflows/local/variant_summary/main'

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



    // Create sampleset files for vcftools pairs analyses
    process WRITE_SAMPLESET {
        tag "$meta.id"
        label 'process_single'

        input:
        tuple val(meta), val(samples)

        output:
        tuple val(meta), val(samples), path("*.txt"),   emit: txt

        script:
        """
        echo -e '${samples.join("\n")}' > ${meta.id}.txt
        """
    }

    ch_sample_sets = WRITE_SAMPLESET(ch_sample_sets).map{
        meta, samples, txt ->
        new_meta = [
            id: meta.id,
            subset: meta.subset,
            samples: samples,
            txt: txt
        ]
        new_meta
    }

    // Build indices
    PREPARE_GENOME(
        fasta,
        fasta_fai
    )
    genome_bed = PREPARE_GENOME.out.genome_bed
    // For now: always add genome bed to ch_roi.
    ch_roi = ch_roi.mix(PREPARE_GENOME.out.genome_bed.map{[[id: "genome", mode: "site"], it[1]]}).mix(genome_bed.map{[[id: "genome", mode: "window"], it[1]]})

    // Gather built indices and resources
    dict               = params.fasta                   ? params.dict                       ? Channel.fromPath(params.dict).collect()                  : PREPARE_GENOME.out.dict                  : []
    fasta_fai          = params.fasta                   ? params.fasta_fai                  ? Channel.fromPath(params.fasta_fai).collect()             : PREPARE_GENOME.out.fasta_fai             : []
    genome_txt         = PREPARE_GENOME.out.genome_txt

    // Build intervals if needed
    PREPARE_INTERVALS(fasta_fai)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    intervals_bed_combined      = params.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_combined  // [interval.bed] all intervals in one file
    intervals                   = PREPARE_INTERVALS.out.intervals_bed        // [interval, num_intervals] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi        = PREPARE_INTERVALS.out.intervals_bed_gz_tbi // [interval_bed, tbi, num_intervals] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather

    // Gather used softwares versions
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
    ch_versions = ch_versions.mix(PREPARE_INTERVALS.out.versions)

    // Summarize coverage for all sample sets. Retrieve mean, var, and
    // sd coverage. Transfer values to auto coverage set.
    MOSDEPTH_COVERAGE(
        ch_input,
        ch_sample_sets,
        fasta,
        fasta_fai,
        genome_txt,
        intervals_bed_combined
    )
    ch_versions = ch_versions.mix(MOSDEPTH_COVERAGE.out.versions.first())

    ch_coverage_set_params = ch_coverage_set.map{
        it -> [[it.id], it]
    }.join(MOSDEPTH_COVERAGE.out.stats.map{
        it -> [[it.id], it]
    }, by: 0, remainder: true).map{
        it, meta, stats ->
        new_meta = [
            id: meta.id,
            subset: meta.subset,
            sampleset_id: meta.sampleset_id,
            coverage_id: meta.coverage_id,
            samples: meta.samples,
            min: meta.min,
            max: meta.max,
        ]
        if (meta.coverage_id == "auto") {
            new_meta.min = [stats.mean - 0.5 * stats.sd, 0].max()
            new_meta.max = stats.mean + 0.5 * stats.sd
        }
        new_meta
    }

    // Given the coverage cutoffs generate coverage masks for
    // sampleset - coverage combinations. coverage_set contains all
    // possible coverage settings and should be combined with the
    // bedgraph file for each sample set
    ch_coverage_set_bed = MOSDEPTH_COVERAGE.out.sum_gz_tbi.map{
        [it[0].sampleset_id, it]
    }.cross(
        ch_coverage_set_params.map{[it.sampleset_id, it]}
    ).map{
        it1, it2 ->
        new_meta = it2[1]
        bed = it1[1][1] // FIXME: make access more readable
        [new_meta, bed]
    }
    genome_txt.view()
    CREATE_COVERAGE_MASKS(ch_coverage_set_bed, genome_txt.map{meta, txt -> txt})
    ch_versions = ch_versions.mix(CREATE_COVERAGE_MASKS.out.versions.first())
    coverage_masks = CREATE_COVERAGE_MASKS.out.gz_tbi

    //
    // SUBWORKFLOW: CREATE_SEQUENCE_MASKS: create sequence (i.e.
    // fasta) mask files from input mask and coverages for
    // single-population analyses

    // FIXME: Add sample-specific present/absent mask, e.g., at least
    // 50% samples in a population must have >0 (or other cutoff)
    // coverage.
    CREATE_SEQUENCE_MASKS (
        PREPARE_GENOME.out.genome_bed,
        fasta,
        ch_roi,
        coverage_masks,
    )
    ch_versions = ch_versions.mix(CREATE_SEQUENCE_MASKS.out.versions.first())

    // Expand mask file list to modes and combine with window sizes
    if (params.window_sizes instanceof Integer) {
        window_sizes = Channel.of(params.window_sizes)
    } else {
        window_sizes = Channel.of(params.window_sizes.split(","))
    }
    CREATE_SEQUENCE_MASKS.out.cov_fasta.branch{
        window: it[0].mode == "window"
        site: it[0].mode == "site"
    }.set{masks_branch}
    masks = masks_branch.window.combine(window_sizes)
        .map{
            meta, fasta, window ->
            new_meta = [
                id: meta.id,
                subset: meta.subset,
                mode: meta.mode,
                sampleset_id: meta.sampleset_id,
                coverage_id: meta.coverage_id,
                samples: meta.samples,
                min: meta.min,
                max: meta.max,
                roi: meta.roi,
                window_size: meta.mode == "window" ? window : null
            ]
            [new_meta, fasta]
        }
        .mix(masks_branch.site)
    VARIANT_SUMMARY(
        vcf.map{[[id:"nucleotide_diversity.${vcf.baseName}"], it]},
        masks
    )

    process WRITE_SAMPLESET {
        tag "$meta.id"
        label 'process_single'

        input:
        tuple val(meta), val(samples)

        output:
        tuple val(meta), path("*.txt"),   emit: txt

        script:
        """
        echo -e '${samples.join("\n")}' > ${meta.id}.txt
        """
    }

    // Create sampleset files for vcftools pairs analyses
    WRITE_SAMPLESET(ch_sample_sets)
    // Need specialized function for pairwise analyses. How to combine
    // masks? E.g. pop1 vs pop2 - intersect the coverage masks?
    masks.branch{
        pair: it[0].subset == true
        single: it[0].subset == false
    }.set{masks_paired}
    //masks_paired.pair.view()

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
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Function to get list of [ meta, bam ]
def create_bam_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample

    // add path(s) of the fastq file(s) to the meta map
    def bam_meta = []
    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> bam file does not exist!\n${row.bam}"
    }
    if (!file(row.bai).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> bai file does not exist!\n${row.bai}"
    }
    bam_meta = [ meta, file(row.bam), file(row.bai) ]
    return bam_meta
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
