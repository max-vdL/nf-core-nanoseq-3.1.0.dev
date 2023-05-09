// Author: Maximilian von der Linde
// Email: maximilian.vonderlinde@ccri.at
// Date: April 3, 2023
// License: MIT

// Description:
// This module includes a Nextflow process that runs deeptools bamCoverage on a set of regions of interest and then plots this coverage for every roi.
// The plotting pyton script is included in the singularity container. It was designed to be integrated into the nannoseq nf-core pipeline.


process COVERAGE_ANALYSIS {
	label 'process_medium'

	conda "bioconda::deeptools bioconda::pyBigWig bioconda::matplotlib bioconda::argparse"
	// container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
	//     'https://depot.galaxyproject.org/singularity/bioconductor-bambu:3.0.8--r42hc247a5b_0' :
	//     'quay.io/biocontainers/bioconductor-bambu:3.0.8--r42hc247a5b_0' }"
	container "$NXF_SINGULARITY_CACHEDIR/deeptools_covplotp-v1.0.img"

	input:
	tuple val(meta), path(sizes), val(is_transcripts), path(input), path(index)
	path fasta
	path bed

	output:
	tuple val(meta), path("*.png"), path("*.pdf"), emit: plot
    path "versions.yml"                    , emit: versions

	when:
	task.ext.when == null || task.ext.when

	script:
	"""
	mkdir bw
	# Loop over regions of interest in bed file
	while read -r chrom start end name; do
		# Define output bigWig file names for each strand
		pos_bigwig="bw/${name::-1}_fwd.bw"
		neg_bigwig="bw/${name::-1}_rev.bw"

		# Generate bigWig files for each strand, excluding reverse-complemented reads
		bamCoverage -b $input -o "$pos_bigwig" -r "chr${chrom}:${start}:${end}" --samFlagExclude 16
		bamCoverage -b $input -o "$neg_bigwig" -r "chr${chrom}:${start}:${end}" --samFlagInclude 16
	done < $bed

	python3 plot_bws.py --bigwig_path bw/ --roi_path $bed -o ./

	cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(deeptools --version 2>&1)
    END_VERSIONS
	"""
}
