process SUM_COVERAGE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed"),    emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outbed = "${prefix}.sum.bed"

    """
    awk  -vFS="\t" -v OFS="\t" 'NR > 1 {
        sum=0;
        for (i=4; i<=NF; i++) sum+=\$i;
        print \$1, \$2, \$3, sum
        }' $bed > $outbed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
