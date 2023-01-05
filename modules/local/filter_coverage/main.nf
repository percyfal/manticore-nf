process FILTER_COVERAGE {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), stdin

    output:
    tuple val(meta), stdout       , emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def min = meta.min ?: 0
    def max = meta.max ?: 1e9
    """
    awk  -vFS="\t" -v OFS="\t" '{
        if (\$4 >= $min && \$4 <= $max) {print }
    }' -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
