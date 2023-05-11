// Author: Maximilian von der Linde
// Email: maximilian.vonderlinde@ccri.at
// Date: Mai 5, 2023
// License: MIT

// Description:
// This module includes a Nextflow process that simply converts one or multiple bed.gz files to ungzipped bed file(s).

process GUNZIPMOSDEPTH {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sed:4.7.0' :
        'quay.io/biocontainers/sed:4.7.0' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
	gunzip $bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | grep "sed (GNU sed)" | sed 's/^.*) //')
    END_VERSIONS
    """
}
