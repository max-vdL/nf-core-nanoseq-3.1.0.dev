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

process CLAIRS {
    tag "$meta.id"
    label 'process_medium'

    container "${NXF_SINGULARITY_CACHEDIR}/hkubal-clairs-0.1.4.img"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(input), path(index)
    tuple path(normal_input), path(normal_index)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("$prefix/output.vcf.gz")    ,  emit: vcf
    tuple val(meta), path("$prefix/output.vcf.gz.tbi"),  emit: tbi  ,  optional: true  
    tuple val(meta), path("$prefix/run_clairs.log") 
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir $prefix
    /opt/bin/run_clairs \\
        --tumor_bam_fn $input \\
        --normal_bam_fn $normal_input \\
        --ref_fn $fasta \\
        --threads $task.cpus \\
        --ctg_name ENST00000318560.6 \\
        --output_dir $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clairs: \$(echo \$(/opt/bin/run_clairs --version 2>&1) | sed 's/^.*run_clairs //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
