/*
 * vcf postprocessing
 */

include { MAPVCF                              } from '../../modules/local/mapvcf'
include { VARIANTANNOTATION                   } from '../../modules/local/variantannotation'
include { VEP2TABLE                   } from '../../modules/local/vep2table'

workflow VCF_POSTPROCESS {
    take:
    ch_vcf
    ch_fasta
    ch_fai

    main:
    ch_mapped_annotated_tables = Channel.empty()
	ch_mapped_annotated_xlsx   = Channel.empty()
    ch_versions                = Channel.empty()

	/*
	* Map the currently to the ensembl reference mapped SNPs to hg38
	*/
	MAPVCF( ch_vcf, params.annotation_hub )
	ch_versions = ch_versions.mix(MAPVCF.out.versions)

	/*
	* Annotate with vep
	*/
	VARIANTANNOTATION( MAPVCF.out.vcf, ch_fasta, params.ensembl_vep )
	ch_versions = ch_versions.mix(VARIANTANNOTATION.out.versions)

	/*
	* Put it into tsv/excel format for nicer viewing
	*/
	VEP2TABLE( VARIANTANNOTATION.out.vcf )
	ch_mapped_annotated_tsv = VEP2TABLE.out.tsv
	ch_mapped_annotated_xlsx = VEP2TABLE.out.xlsx
	ch_versions = ch_versions.mix(VEP2TABLE.out.versions)

    emit:
    ch_mapped_annotated_tsv
    ch_mapped_annotated_xlsx
    ch_versions
}