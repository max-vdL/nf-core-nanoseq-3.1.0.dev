// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process VARIANTANNOTATION {
    tag "$meta.id"
    label 'process_single'

    conda "ensembl-vep=103.1"
    container "${NXF_SINGULARITY_CACHEDIR}/ensembl-vep-v103.1.img"

    input:
    tuple val(meta), path(vcf)
    path fasta
	path ensembl_vep 

    output:
    tuple val(meta), path("*.annot.vcf")          , emit: vcf
    tuple val(meta), path("*.annot.vepstats.html"), emit: stats
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Annotating vcf file: $vcf"

    input_file=${vcf}
    output_file=\$(basename $vcf .vcf).filt.gmap.annot.vcf
    stats_file=\$(basename $vcf .vcf).filt.gmap.annot.vepstats.html

    vep --fork $task.cpus \
        --force_overwrite \
        --fasta $fasta \
        --input_file $vcf \
        --format vcf \
        --output_file \$output_file \
        --vcf \
        --stats_file \$stats_file \
        --species homo_sapiens \
        --assembly GRCh38 \
        --offline \
        --cache \
        --cache_version 103 \
        --everything \
        --dir_cache $ensembl_vep


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembl-vep: '103.1'
    END_VERSIONS
    """
}
