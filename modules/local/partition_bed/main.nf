process PARTITION_BED {
    tag "$intervals"

    conda (params.enable_conda ? "conda-forge::python=3.9.15 pandas=1.4.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'quay.io/biocontainers/pandas:1.4.3' }"

    input:
    path(intervals)

    output:
    path("*.bed")         , emit: bed
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    partition_bed.py \\
        $intervals \\
        --npartitions ${params.num_intervals}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
