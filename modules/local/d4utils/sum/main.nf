process D4UTILS_SUM {
    tag "$meta.id"
    label 'process_single'

    conda "${projectDir}/modules/local/d4utils/d4utils.yml"
    // FIXME: no container available for d4tools
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'quay.io/biocontainers/YOUR-TOOL-HERE' }"


    input:
    tuple val(meta), path(d4)

    output:
    tuple val(meta), path("*.d4"), emit: d4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //def input_files = label.join(' ')

    // FIXME: add intervals argument and output to file name
    // corresponding to interval name. Then we collect the intermediate results
    """
    d4utils.py \\
        sum \\
        $args \\
        $d4 \\
        ${prefix}.d4
    """
}
