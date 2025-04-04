# ![nf-core/nanoseq](docs/images/nf-core-nanoseq_logo_light.png#gh-light-mode-only) ![nf-core/nanoseq](docs/images/nf-core-nanoseq_logo_dark.png#gh-dark-mode-only)

[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.3697959-1073c8)](https://doi.org/10.5281/zenodo.3697959)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)


## Introduction

**CCRI/nfcore/nanoseq-3.1.0.dev** is a bioinformatics analysis pipeline for Nanopore DNA/RNA sequencing data. It is an adaptation of the nf-core nanoseq pipeline (version 3.1.0) with additional optional functionalities. This adaptation retains all functionalities from the original nf-core pipeline. Any parameter files generated with the original pipeline should still work in this adaptation.


### Changes Made
The following optional additions have been made to the pipeline:

- Multiple file input: A sample can now consist of multiple files of the same file type (e.g., fastq.gz) using wildcard paths.
- Multiple variant callers can now be specified and are run in paralell.
- Structural variant callers: The pipeline includes the following structural variant callers:
  - [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)
  - [Dysgu](https://github.com/kcleal/dysgu)
  - [DeBreak](https://github.com/Maggi-Chen/DeBreak)
- Read depth analysis: The pipeline includes [Mosdepth](https://github.com/brentp/mosdepth) for read depth analysis in regions of interest. A custom module was created. 
  - Connected to this is the custom module called gunzipmosdepth which unzips the mosdepth output. A custom module was created.
- Custom read depth plotting module: A custom read depth plotting module called `coverageanalysis` has been added. Please note that this script functions only with its associated Singularity/Docker image (REF) and cannot be run in Conda mode. To specify regions of interest, provide a bed file.
- SNP/INDEL caller: The pipeline includes the following SNP/INDEL callers:
  - [Clair3](https://github.com/HKU-BAL/Clair3) for ML based germline mutation calling, used nf-core module.
  - [ClairS](https://github.com/HKU-BAL/ClairS) for ML based somatic mutation calling, used nf-core module.
  - [Medaka-v1.8.0](https://github.com/nanoporetech/medaka/releases/tag/v1.8.0) for gvcf output which the implemented Medaka version can not provide, created custom module.
  - [Deepvariant-v1.5.0](https://github.com/google/deepvariant/releases/tag/v1.5.0)
- SNP variant calling postprocessing: The pipeline incorporates some functions from the [nanopanel2](https://github.com/popitsch/nanopanel2) pipeline for SNP variant calling postprocessing. The postprocessing steps include:
  - Mapping variants to a UCSC chromosome (e.g., hg38) if alignment was performed with an Ensembl transcript as the reference. The nf-core module mapvcf was used.
  - Annotation with [VEP](https://github.com/Ensembl/ensembl-vep). Used the nf-core module variantannotation.
  - Transferral of variants into TSV and Excel formats with custom module vep2table.
- UMI handling with [umitools](https://github.com/CGATOxford/UMI-tools/releases/tag/1.1.4) dedup and extract: 

## nf-core/nanoseq Documentation
For information on functionality, documentation and anything else taken over from nf-core nanoseq, please visit their [Github](https://github.com/nf-core/nanoseq/tree/3.1.0)