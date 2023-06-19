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

process CLAIR3 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::clair3=1.0.2"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/clair3:1.0.2--py39hb9dc472_0':
    //     'quay.io/biocontainers/clair3:1.0.2--py39hb9dc472_0' }"
    container "${NXF_SINGULARITY_CACHEDIR}/clair3_latest.img"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(input), path(index)
    path(fasta)
    path(fai)


    output:
    tuple val(meta), path("$prefix/merge_output.vcf.gz")    ,  emit: vcf
    tuple val(meta), path("$prefix/merge_output.gvcf.gz")   ,  emit: gvcf  ,  optional: true
    tuple val(meta), path("$prefix/merge_output.*vcf.gz.tbi"),  emit: tbi  ,  optional: true  
    path "versions.yml"                  , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"

    """
    run_clair3.sh \\
        --bam_fn=$input \\
        --ref_fn=$fasta \\
        --threads=$task.cpus \\
        --platform="ont" \\
        --model_path="/opt/models/r941_prom_sup_g5014" \\
        --output=$prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(run_clair3.sh --version 2>&1) | sed 's/^.*Clair3 v//; s/Using.*\$//' ))
    END_VERSIONS
    """
}
