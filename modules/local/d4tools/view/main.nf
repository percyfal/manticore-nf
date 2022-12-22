process D4TOOLS_VIEW {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::d4tools=0.3.8"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(d4)


    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    d4tools \\
        view \\
        $args \\
        $d4 \\
        > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        d4tools: \$(echo \$(d4tools view --version 2>&1) | sed 's/d4-view - View a D4 file //;' )
    END_VERSIONS
    """
}
