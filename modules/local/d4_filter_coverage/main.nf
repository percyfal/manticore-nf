process D4_FILTER_COVERAGE {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::gawk=5.1.0 bioconda::d4tools==0.3.8" : null)
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
    //    'quay.io/biocontainers/gawk:5.1.0' }"


    input:
    tuple val(meta), path(d4)

    output:
    tuple val(meta), path("*.bed"),    emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed = "${prefix}.sum.bed"
    def samples = "\\(${meta.samples.join('\\|')}\\)"
    """
    columns=\$(d4tools ls-track $d4 | grep -n "$samples" | sed "s/:[^ ]*//g" | awk '{print \$1 + 3 }' | tr "\\n" "," | sed "s/,\$//")
    d4tools \\
        view \\
        $d4 | cut -f "1,2,3,\$columns" | awk  -vFS="\t" -v OFS="\t" '{
        sum=0;
        for (i=4; i<=NF; i++) sum+=\$i;
        print \$1, \$2, \$3, sum
        }' > $bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
