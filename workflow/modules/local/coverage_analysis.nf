process COVERAGE_ANALYSIS {
    label 'process_medium'

    // conda "conda-forge::r-base=4.0.3 bioconda::bioconductor-bambu=3.0.8 bioconda::bioconductor-bsgenome=1.66.0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/bioconductor-bambu:3.0.8--r42hc247a5b_0' :
    //     'quay.io/biocontainers/bioconductor-bambu:3.0.8--r42hc247a5b_0' }"
	container "$NXF_SINGULARITY_CACHEDIR/deeptools_covplotp-v1.0.img"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(input), path(index)
    path bed
	path fasta

    output:

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
    """
}
