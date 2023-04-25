// Custom module, NOT TESTED! Was discontinued in favor of nf-core module mosdepth

process MOSDEPTH {
    label 'process_medium'

    conda "bioconda::mosdepth=0.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.3--h37c5b7d_2' :
        'quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2' }"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(input), path(index)
    path bed


    output:
    path "*dist.txt"   			   				, emit: dist
    path "*summary.txt"   			   			, emit: summary
    tuple path("*per-base*"), path("*regions*") , emit: inputs
    path "versions.yml"                    		, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mosdepth \
        -b $bed \
        ${meta} \
        ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles: \$(sniffles --help 2>&1 | head -n 1 //')
    END_VERSIONS
    """
}

