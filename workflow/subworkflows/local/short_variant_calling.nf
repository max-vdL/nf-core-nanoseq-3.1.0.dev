/*
 * Short variant calling test
 */

include { MEDAKA_VARIANT                        } from '../../modules/local/medaka_variant_v1.8.0'
include { TABIX_BGZIP as MEDAKA_BGZIP_VCF       } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as MEDAKA_BGZIP_GVCF       } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as MEDAKA_TABIX_VCF       } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as MEDAKA_TABIX_GVCF       } from '../../modules/nf-core/tabix/tabix/main'
include { DEEPVARIANT                           } from '../../modules/local/deepvariant'
include { TABIX_TABIX as DEEPVARIANT_TABIX_VCF  } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as DEEPVARIANT_TABIX_GVCF } from '../../modules/nf-core/tabix/tabix/main'
include { PEPPER_MARGIN_DEEPVARIANT             } from '../../modules/local/pepper_margin_deepvariant'
include { CLAIR3                                } from '../../modules/local/clair3'
include { CLAIRS                                } from '../../modules/local/clairs'
include { TABIX_TABIX as CLAIR3_TABIX_VCF  } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as CLAIR3_TABIX_GVCF } from '../../modules/nf-core/tabix/tabix/main'

workflow SHORT_VARIANT_CALLING {

    take:
    ch_view_sortbam
    ch_fasta
    ch_fai

    main:
    ch_short_calls_vcf              = Channel.empty()
    ch_short_calls_vcf_tbi          = Channel.empty()
    ch_short_calls_gvcf             = Channel.empty()
    ch_short_calls_gvcf_tbi         = Channel.empty()
    ch_versions                     = Channel.empty()

    /*
     * Call short variants
     */
    if (params.variant_caller.split(',').contains('medaka')) {

        /*
         * Call short variants with medaka
         */
        MEDAKA_VARIANT( ch_view_sortbam, ch_fasta )
        ch_versions = ch_versions.mix(medaka_version = MEDAKA_VARIANT.out.versions)

        /*
         * Zip medaka vcf
         */
        MEDAKA_BGZIP_VCF( MEDAKA_VARIANT.out.vcf )
        ch_short_calls_vcf  = MEDAKA_BGZIP_VCF.out.output
        ch_versions = ch_versions.mix(bgzip_version = MEDAKA_BGZIP_VCF.out.versions)
        /*
         * Zip medaka g.vcf
         */
        // MEDAKA_BGZIP_GVCF( MEDAKA_VARIANT.out.gvcf )
        // ch_short_calls_gvcf  = MEDAKA_BGZIP_GVCF.out.output

        /*
         * Index medaka vcf.gz
         */
        MEDAKA_TABIX_VCF( ch_short_calls_vcf )
        ch_short_calls_vcf_tbi  = MEDAKA_TABIX_VCF.out.tbi
        ch_versions = ch_versions.mix(tabix_version = MEDAKA_TABIX_VCF.out.versions)
        /*
         * Index medaka g.vcf.gz
         */
        // MEDAKA_TABIX_GVCF( ch_short_calls_gvcf )
        // ch_short_calls_gvcf_tbi  = MEDAKA_TABIX_GVCF.out.tbi

    }
	if (params.variant_caller.split(',').contains('deepvariant')) {

        /*
        * Call variants with deepvariant
        */
        DEEPVARIANT( ch_view_sortbam, ch_fasta, ch_fai )
        ch_short_calls_vcf  = DEEPVARIANT.out.vcf
        ch_short_calls_gvcf = DEEPVARIANT.out.gvcf
        ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)

        /*
         * Index deepvariant vcf.gz
         */
        DEEPVARIANT_TABIX_VCF( ch_short_calls_vcf )
        ch_short_calls_vcf_tbi  = DEEPVARIANT_TABIX_VCF.out.tbi
        ch_versions = ch_versions.mix(DEEPVARIANT_TABIX_VCF.out.versions)

        /*
         * Index deepvariant g.vcf.gz
         */
        DEEPVARIANT_TABIX_GVCF( ch_short_calls_gvcf )
        ch_short_calls_gvcf_tbi  = DEEPVARIANT_TABIX_GVCF.out.tbi
        ch_versions = ch_versions.mix(DEEPVARIANT_TABIX_VCF.out.versions)

    }
	if (params.variant_caller.split(',').contains('clair3')) {

        /*
        * Call variants with clair3
        */
        CLAIR3( ch_view_sortbam, ch_fasta, ch_fai )
        ch_short_calls_vcf  = CLAIR3.out.vcf
        // ch_short_calls_gvcf = CLAIR3.out.gvcf
        ch_versions = ch_versions.mix(CLAIR3.out.versions)

        /*
         * Index clair3 vcf.gz
         */
        // CLAIR3_TABIX_VCF( ch_short_calls_vcf )
        // ch_short_calls_vcf_tbi  = CLAIR3_TABIX_VCF.out.tbi
        // ch_versions = ch_versions.mix(CLAIR3_TABIX_VCF.out.versions)

        // /*
        //  * Index clair3 g.vcf.gz
        //  */
        // CLAIR3_TABIX_GVCF( ch_short_calls_gvcf )
        // ch_short_calls_gvcf_tbi  = CLAIR3_TABIX_GVCF.out.tbi
        // ch_versions = ch_versions.mix(CLAIR3_TABIX_GVCF.out.versions)

    } 
	if (params.variant_caller.split(',').contains('clairs')) {
		channel.fromPath(params.clairs_normbam)
			.map { file -> tuple(file, file + '.bai') }
			.set { ch_clairs_normbam }        
		CLAIRS( ch_view_sortbam, ch_clairs_normbam, ch_fasta, ch_fai )

        ch_short_calls_vcf  = CLAIRS.out.vcf
        ch_versions = ch_versions.mix(CLAIRS.out.versions)
	}
	if (params.variant_caller.split(',').contains('pepper_margin_deepvariant')) {

        /*
         * Call variants with pepper_margin_deepvariant (automatic zip + index, docker + singularity only)
         */
        PEPPER_MARGIN_DEEPVARIANT( ch_view_sortbam, ch_fasta, ch_fai )
        ch_short_calls_vcf = PEPPER_MARGIN_DEEPVARIANT.out.vcf
        ch_short_calls_vcf_tbi = PEPPER_MARGIN_DEEPVARIANT.out.tbi
        ch_versions = ch_versions.mix(PEPPER_MARGIN_DEEPVARIANT.out.versions)
    }

    emit:
    ch_short_calls_vcf
    ch_short_calls_vcf_tbi
    ch_short_calls_gvcf
    ch_short_calls_gvcf_tbi
    ch_versions
}
