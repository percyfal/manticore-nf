process BUILD_INTERVALS {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(fasta_fai)

    output:
    tuple val(meta), path("${fasta_fai.baseName}.bed") , emit: bed
    tuple val(meta), path("${fasta_fai.baseName}.txt") , emit: txt
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fasta_fai} > ${fasta_fai.baseName}.bed
    awk -v FS='\t' -v OFS='\t' '{ print \$1, \$2 }' ${fasta_fai} > ${fasta_fai.baseName}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
