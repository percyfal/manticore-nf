process D4TOOLS_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::d4tools=0.3.8"
    // FIXME: no container available for d4tools
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(d4), val(label)

    output:
    tuple val(meta), path("*.d4"), emit: d4
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_files = label.join(' ')

    // FIXME: add intervals argument and output to file name
    // corresponding to interval name. Then we collect the intermediate results
    """
    d4tools \\
        merge \\
        $args \\
        $input_files \\
        ${prefix}.d4

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        d4tools: \$(echo \$(d4tools merge --version 2>&1) | sed 's/d4-merge - Merge multiple D4 files //;' )
    END_VERSIONS
    """
}
