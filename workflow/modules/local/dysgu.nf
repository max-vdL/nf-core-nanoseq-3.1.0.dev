// Author: Maximilian von der Linde
// Email: maximilian.vonderlinde@ccri.at
// Date: April 3, 2023
// License: MIT

// Description:
// This module includes a Nextflow process that runs the SV caller dysgu (https://github.com/kcleal/dysgu) in nanopore mode. 
// It was designed to be integrated into the nannoseq nf-core pipeline.

process DYSGU {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::dysgu=1.3.16"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dysgu%3A1.3.16--py310h0de0465_0' :
        'quay.io/repository/biocontainers/dysgu:1.3.16--py310h0de0465_0' }"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(input), path(index)
    path(fasta)


    output:
    tuple val(meta), path("*_dysgu.vcf")   , emit: sv_calls
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    dysgu call \\
        --mode nanopore \\
        --pl nanopore \\
        -p $task.cpus \\
        -o ${meta.id}_dysgu.vcf \\
        ${args} \\
        ${fasta} \\
        \$PWD/tmp \\
        ${input} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dysgu: \$(dysgu --version 2>&1 | grep version | sed 's/, version//')
    END_VERSIONS
    """
}

