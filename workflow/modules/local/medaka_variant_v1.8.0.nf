process MEDAKA_VARIANT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::medaka=1.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.8.0--py310hb7fe8e6_0' :
        'quay.io/biocontainers/medaka:1.8.0--py310hb7fe8e6_0' }"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(input), path(index)
    path(fasta)

    output:
    tuple val(meta), path ("$output_vcf"), emit: vcf // vcf files
    tuple val(meta), path ("$output_gvcf"), emit: gvcf // gvcf files
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    output_dir = "${meta.id}"
    output_vcf = output_dir+"/*.vcf"
    output_gvcf = output_dir+"/*.g.vcf"
    """

    medaka variant \\
        $fasta \\
        $input \\
        $output_dir \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
