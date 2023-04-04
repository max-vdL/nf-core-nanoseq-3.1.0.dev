// Author: Maximilian von der Linde
// Email: maximilian.vonderlinde@ccri.at
// Date: April 3, 2023
// License: MIT

// Description:
// This module includes a Nextflow process that runs the SV caller DeBreak (https://github.com/Maggi-Chen/DeBreak) in nanopore mode. 
// It was designed to be integrated into the nannoseq nf-core pipeline.

process DEBREAK {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::debreak=1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/debreak:1.3--h9ee0642_0' :
        'quay.io/repository/biocontainers/debreak:1.3--h9ee0642_0' }"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(input), path(index)
    path(fasta)


    output:
    tuple val(meta), path("*_debreak.vcf")   , emit: sv_calls
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    debreak \
        -t $task.cpus \
		--bam ${input} \
        -o ${meta.id}_debreak.vcf \
        -r ${fasta} \
		--rescue_large_ins --rescue_dup --poa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        debreak: \$(debreak --version 2>&1 | grep version | sed 's/version//')
    END_VERSIONS
    """
}

