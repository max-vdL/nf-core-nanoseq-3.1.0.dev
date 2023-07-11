
process VEP2TABLE {
    tag "$meta.id"
    label 'process_single'

    container "${NXF_SINGULARITY_CACHEDIR}/annot_table_v1.1.img"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.xlsx"                 , emit: xlsx
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vepClaire3Vcf2table.r --file=$vcf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep2table: '1.0')
    END_VERSIONS
    """
}
