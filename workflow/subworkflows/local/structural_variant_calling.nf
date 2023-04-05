/*
 * Structural variant calling
 */

include { SNIFFLES                              } from '../../modules/local/sniffles'
include { SNIFFLES2                             } from '../../modules/nf-core/sniffles'
include { BCFTOOLS_SORT as SNIFFLES_SORT_VCF    } from '../../modules/nf-core/bcftools/sort/main'
include { TABIX_BGZIP as SNIFFLES_BGZIP_VCF     } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as SNIFFLES_TABIX_VCF     } from '../../modules/nf-core/tabix/tabix/main'
include { CUTESV                                } from '../../modules/local/cutesv'
include { BCFTOOLS_SORT as CUTESV_SORT_VCF      } from '../../modules/nf-core/bcftools/sort/main'
include { TABIX_BGZIP as CUTESV_BGZIP_VCF       } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as CUTESV_TABIX_VCF       } from '../../modules/nf-core/tabix/tabix/main'
include { DYSGU                                } from '../../modules/local/dysgu'
include { BCFTOOLS_SORT as DYSGU_SORT_VCF      } from '../../modules/nf-core/bcftools/sort/main'
include { TABIX_BGZIP as DYSGU_BGZIP_VCF       } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as DYSGU_TABIX_VCF       } from '../../modules/nf-core/tabix/tabix/main'
include { DEBREAK                                } from '../../modules/local/debreak'
include { BCFTOOLS_SORT as DEBREAK_SORT_VCF      } from '../../modules/nf-core/bcftools/sort/main'
include { TABIX_BGZIP as DEBREAK_BGZIP_VCF       } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as DEBREAK_TABIX_VCF       } from '../../modules/nf-core/tabix/tabix/main'


workflow STRUCTURAL_VARIANT_CALLING {

    take:
    ch_view_sortbam
    ch_fasta
    ch_fai

    main:
    ch_sv_calls_vcf     = Channel.empty()
    ch_sv_calls_vcf_tbi = Channel.empty()

    ch_versions         = Channel.empty()

    /*
     * Call structural variants with sniffles
     */
    if (params.structural_variant_caller == 'sniffles') {

        /*
         * Call structural variants with sniffles
         */
        SNIFFLES( ch_view_sortbam )
        ch_versions = ch_versions.mix(SNIFFLES.out.versions)

        /*
         * Sort structural variants with bcftools
         */
        SNIFFLES_SORT_VCF( SNIFFLES.out.sv_calls )
        ch_sv_calls_vcf = SNIFFLES_SORT_VCF.out.vcf
        ch_versions = ch_versions.mix(SNIFFLES_SORT_VCF.out.versions)

        /*
         * Index sniffles vcf.gz
         */
        SNIFFLES_TABIX_VCF( ch_sv_calls_vcf )
        ch_sv_calls_tbi  = SNIFFLES_TABIX_VCF.out.tbi
        ch_versions = ch_versions.mix(SNIFFLES_TABIX_VCF.out.versions)

    } else if (params.structural_variant_caller == 'sniffles2') {

        /*
         * Call structural variants with sniffles2
         */
        SNIFFLES2( ch_view_sortbam, ch_fasta )
        ch_versions = ch_versions.mix(SNIFFLES2.out.versions)

        /*
         * Sort structural variants with bcftools
         */
        SNIFFLES_SORT_VCF( SNIFFLES2.out.sv_calls )
        ch_sv_calls_vcf = SNIFFLES_SORT_VCF.out.vcf
        ch_versions = ch_versions.mix(SNIFFLES_SORT_VCF.out.versions)

        /*
         * Index sniffles2 vcf.gz
         */
        SNIFFLES_TABIX_VCF( ch_sv_calls_vcf )
        ch_sv_calls_tbi  = SNIFFLES_TABIX_VCF.out.tbi
        ch_versions = ch_versions.mix(SNIFFLES_TABIX_VCF.out.versions)

    } else if (params.structural_variant_caller == 'cutesv') {

        /*
        * Call structural variants with cutesv
        */
        CUTESV( ch_view_sortbam, ch_fasta )
        ch_versions = ch_versions.mix(CUTESV.out.versions)

        /*
         * Sort structural variants with bcftools
         */
        CUTESV_SORT_VCF( CUTESV.out.sv_calls )
        ch_sv_calls_vcf = CUTESV_SORT_VCF.out.vcf
        ch_versions = ch_versions.mix(CUTESV_SORT_VCF.out.versions)

        /*
         * Zip cutesv vcf.gz
         */
        CUTESV_TABIX_VCF( ch_sv_calls_vcf )
        ch_sv_calls_tbi  = CUTESV_TABIX_VCF.out.tbi
        ch_versions = ch_versions.mix(CUTESV_TABIX_VCF.out.versions)

    } else if (params.structural_variant_caller == 'dysgu') {

        /*
        * Call structural variants with dysgu
        */
        DYSGU( ch_view_sortbam, ch_fasta )
        ch_versions = ch_versions.mix(DYSGU.out.versions)

        /*
         * Sort structural variants with bcftools
         */
        DYSGU_SORT_VCF( DYSGU.out.sv_calls )
        ch_sv_calls_vcf = DYSGU_SORT_VCF.out.vcf
        ch_versions = ch_versions.mix(DYSGU_SORT_VCF.out.versions)

        /*
         * Zip dysgu vcf.gz
         */
        DYSGU_TABIX_VCF( ch_sv_calls_vcf )
        ch_sv_calls_tbi  = DYSGU_TABIX_VCF.out.tbi
        ch_versions = ch_versions.mix(DYSGU_TABIX_VCF.out.versions)

    } else if (params.structural_variant_caller == 'debreak') {

        /*
        * Call structural variants with debreak
        */
        DEBREAK( ch_view_sortbam, ch_fasta )
        ch_versions = ch_versions.mix(DEBREAK.out.versions)

        /*
         * Sort structural variants with bcftools
         */
        DEBREAK_SORT_VCF( DEBREAK.out.sv_calls )
        ch_sv_calls_vcf = DEBREAK_SORT_VCF.out.vcf
        ch_versions = ch_versions.mix(DEBREAK_SORT_VCF.out.versions)

        /*
         * Zip debreak vcf.gz
         */
        DEBREAK_TABIX_VCF( ch_sv_calls_vcf )
        ch_sv_calls_tbi  = DEBREAK_TABIX_VCF.out.tbi
        ch_versions = ch_versions.mix(DEBREAK_TABIX_VCF.out.versions)
    }

    emit:
    ch_sv_calls_vcf
    ch_sv_calls_vcf_tbi
    ch_versions
}
