process COVERAGE_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.stats.txt")                       ,      emit: txt
    tuple val(meta), stdout,                              emit: stats
    path "versions.yml"              ,                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stats = "${prefix}.stats.txt"

    """
    outstats=\$(awk  -vFS="\t" -v OFS="\t" '
        BEGIN  {sum=0; sq=0; total=0;}
        { sum+=\$4 * (\$3 - \$2); sq+=\$4*\$4 * (\$3 - \$2); total+=(\$3 - \$2);}
        END {
            var=sq/(total) - (sum/total)**2;
            printf("%i\\t%f\\t%f\\t%f\\n", total, sum/total, var, sqrt(var));
        }' $bed
    )
    echo -e "records\\tmean\\tvar\\tsd" >$stats
    echo "\$outstats" >> $stats
    echo -n "\$outstats"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
