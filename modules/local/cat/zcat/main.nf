process CAT_ZCAT {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::pigz=2.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'quay.io/biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(infile)

    output:
    tuple val(meta), stdout         , emit: txt
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    // | input     | output     | command1 | command2 |
    // |-----------|------------|----------|----------|
    // | gzipped   | unzipped   | zcat      |          |
    // | ungzipped | ungzipped  | cat      |          |

    // Use input file ending as default
    prefix   = task.ext.prefix ?: "${meta.id}"
    in_zip   = infile.toString().endsWith('.gz')
    command = in_zip ? 'zcat' : 'cat'
    """
    $command \\
        $args \\
        ${infile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
