// Author: Maximilian von der Linde
// Email: maximilian.vonderlinde@ccri.at
// Date: Mai 5, 2023
// License: MIT

// Description:
// This module includes a Nextflow process that runs deeptools bamCoverage on a set of regions of interest and then plots this coverage for every roi.
// The plotting pyton script is included in the singularity container. It was designed to be integrated into the nannoseq nf-core pipeline.

process COVERAGEANALYSIS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::deeptools=3.5.1 bioconda::pybigwig=0.3.18 conda-forge::matplotlib=3.3.4 bioconda::argparse=1.4.0"
	container "${NXF_SINGULARITY_CACHEDIR}/deeptools_covplot-v1.1.img"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(input), path(index)
    tuple val(meta2), path(fasta)
	path bed

    output:
	tuple val(meta), path("*.png"), path("*.pdf"),  emit: plot
	tuple val(meta), path("bw/")  ,                 emit: bws
    path "versions.yml"           ,                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
	mkdir bw
	# Loop over regions of interest in bed file
	while read -r chrom start end name; do
		# remove newline from last column
		name=`echo "\$name" | sed 's/.\$//'`
		# Define output bigWig file names for each strand
		pos_bigwig="bw/\${name}_fwd.bw"
		neg_bigwig="bw/\${name}_rev.bw"

		# Generate bigWig files for each strand, excluding reverse-complemented reads
		bamCoverage -b $input -o "\$pos_bigwig" -r "chr\${chrom}:\${start}:\${end}" --samFlagExclude 16
		bamCoverage -b $input -o "\$neg_bigwig" -r "chr\${chrom}:\${start}:\${end}" --samFlagInclude 16
	done < $bed

	python3 /scripts/plot_bws.py --bigwig_path bw/ --roi_path $bed -o "./${meta.id}_coverage_plot"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(deeptools --version 2>&1 | sed 's/deeptools //')
    END_VERSIONS
    """
}
