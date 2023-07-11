#!/usr/bin/env Rscript

# converting vep annotated vcf into csv
# re-write in python; compare to Niko's and/or Matthias script (if they are using VEP) or to annovar!
# extracting important fields
# [ ] - considering only filtered and annotated for now!
# [ ] - check results in IGV
# [ ] - add option to specify csv/tsv
# [ ] - add conflict function to highlight any potential conflicts to be resolved by import 
# [ ] - fix issue with Error: ALT length != 1 - is the vcf file normalized?
# [ ] - fix issue - raid issue?
# [ ] - add unit tests and check when run on command line that it produces desired output!
#In open.connection(file, "wb") :
#  cannot open file '/raid/nanopanel2/np2_results/Human_GRCh38_v103/no_truth/flongle20_GRCh38/flongle20_GRCh38.S06_genomic_filt_annot.mm2.tsv': No such file or directory
#Execution halted

# importing key functions using import ----
import::from(.from = optparse, make_option, OptionParser, parse_args, print_help)
import::from(.from = cli, cli_rule, cli_progress_bar, cli_alert_info, cli_alert_warning, cli_progress_step, cat_line, cli_alert_success, cli_progress_done)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="vcf file in genomic coordinates and vep annotated - *.filt.gmap.annot.vcf", metavar="character")) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# checking input files ----
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Input vcf file needs to be supplied!", call.=FALSE)
}

# importing additional functions ----
# not needed to import if the above check fails
import::from(.from = VariantAnnotation, readVcf, isSNV, isIndel, header, info, geno, meta) # genotypeToSnpMatrix
import::from(ensemblVEP, parseCSQToGRanges) 
import::from(.from = IRanges, CharacterList)
import::from(.from = S4Vectors, DataFrame, unstrsplit) 
import::from(.from = BiocGenerics, as.data.frame)
import::from(.from = magrittr, "%>%")
import::from(.from = dplyr, select, rename, left_join, all_of) # filter
import::from(.from = SummarizedExperiment, rowRanges)
import::from(.from = readr, write_tsv)
import::from(.from = openxlsx, write.xlsx)
import::from(.from = stringr, str_split)

# to help parsing the CSQ field
# https://bioconductor.org/packages/release/bioc/html/ensemblVEP.html
#BiocManager::install("ensemblVEP")

cli_rule(left = "vepVcf2table")
cli_progress_bar("Converting annotated vcf to table", type="tasks")

#suppressMessages(library("VariantAnnotation"))  
#library("tidyverse")

# define input_vcf ----
input_vcf <- file.path(opt$file)

cli_progress_step("Loading vcf file.")
if (!file.exists(input_vcf)){
  stop("vcf file does not exists!", call.=FALSE)
} else {
  cli_alert_info("Processing: {input_vcf}")
}

# define output vcf
output_basename <- gsub(pattern = ".vcf$", replacement = "", x=input_vcf)

# 220509 Somehow readVcf has a problem with the vep processed files from rule annotate.
# So I instead use a custom txt file handler for empty vcf files.
tmp_vcf <- readLines(input_vcf)
tmp_vcf <- tmp_vcf[!grepl("#", tmp_vcf)]
if(length(tmp_vcf)==0){
  df <- data.frame(matrix(ncol = 34, nrow = 0))
  colnames(df) <- c("ID", "SYMBOL", "REF", "ALT", "VAR_RAW", "VAF_CORR", "CNT_RAW", "CNT_CORR", "Protein_position", "Amino_acids", "Codons",
                  "HGVSp", "Existing_variation",
                  "Consequence", "IMPACT", "SIFT", "PolyPhen",
                  "HGVSc", "CHROM", "POS", "MANE_SELECT","Gene", "EXON", "cDNA_position",	"CDS_position",	 
                  "QUAL", "BASE_QUAL", "FILTER", "TYPE", "CTX",
                  "SB_P", "AA_SKEW", "HP_FRAC",	"HP_LEN")
  output_tsv <- paste0(output_basename, ".tsv")
  output_excel <- paste0(output_basename, ".xlsx")
  cli_progress_step("Saving empty table as tsv and excel file.")
  cli_alert_info("tsv file: {output_tsv}")
  cli_alert_info("excel file: {output_excel}")

  readr::write_tsv(df, file = output_tsv)
  openxlsx::write.xlsx(df, file = output_excel, asTable = TRUE, overwrite = TRUE)

  cli_alert_success("No need to process VCF.")
  cli_progress_done()
} else {
	# using VariantAnnotation ----
	#message("loading vcf file...")
	vcf <- VariantAnnotation::readVcf(input_vcf, genome = "GRCh38")
	cli_alert_info("Variants: {nrow(vcf)} SNVs: {sum(isSNV(vcf))} INDELs: {sum(isIndel(vcf))}")

	# extracting key fields ----
	# [ ] - positions should match, but double check!
	vcf_header <- header(vcf)
	vcf_gr <- rowRanges(vcf)
	vcf_info <- info(vcf)
	vcf_geno <- geno(vcf)
	# TO-DO
	# combine rowRanges(vcf), vcf_info_df (-CSQ) and CSQ splitted!; remove empty fields from splitted! (potentially keep?!)
	cli_progress_step("Processing GRanges.")
	vcf_gr$ID <- names(vcf_gr) # preserving vcf _name original name
	# just in case to convert REF = DNAStringSet, ALT = DNAStringSetList to character
	# https://support.bioconductor.org/p/66874/
	# ?? https://support.bioconductor.org/p/83185/

	if (any(lapply(vcf_gr$ALT, length) != 1)){
		# checking if any of the ALT contains other variants
		# [ ] - find a test example!
		stop("ALT length != 1 - is the vcf file normalized?", call.=FALSE)
	} else {
		vcf_gr$REF <- as.character(vcf_gr$REF)
		vcf_gr$ALT <- unstrsplit(CharacterList(vcf_gr$ALT), sep = ",")  # as.character(unlist(vcf_gr$ALT)), but may not work with multiple SNVs!
	}

	#cli::cat_boxx(paste(names(vcf_gr), collapse = "_"))
	# additional option to convert gr to data.frame: https://stackoverflow.com/questions/59370461/granges-as-column-in-basedata-frame
	# default base::as.data.frame does not convert gr names to rownames!
	# BiocGenerics::as.data.frame
	vcf_ranges_df <- as.data.frame(vcf_gr) %>%
		rename(CHROM = seqnames,POS = start) %>%
		select(-end, -width, -strand, -paramRangeID)

	# extracting vep info, separating and adding column names from
	# TODO: find easier/better way!
	# vcf@metadata$header@header@listData$INFO
	# ! Limit VEP only to some values!
	cli_progress_step("Processing Info and CSQ fields.")
	#vcf_header_info <- info(vcf_header)

	# pre-select columns
	# rename according to Niko's naming!
	# e.g. AF = vaf_raw
	# PASS_CALLS FILTER_CALLS ALL_CALLS - applicable to consensus calls; not here
	# BASE_QUAL  = avg_base_qual
	# STATUS - truth status (remove as not applicable?) - but keep if it is        
	#   write.table(as.data.frame(vcf_info), file='/home/max_vdl/bioinf_isilon/core_bioinformatics_unit/Internal/max_vdl/labdia/nanopanel2/flongle_39/tmp.tsv')
	if (length(vcf_info) >= 19) {
		vcf_info_df <- as.data.frame(vcf_info) %>%
			rename(VAR_RAW = AF) %>%
			select(-PASS_CALLS, -FILTER_CALLS, -ALL_CALLS, -CSQ)

		keep_columns <- c("ID", "SYMBOL", "REF", "ALT", "VAR_RAW", "VAF_CORR", "CNT_RAW", "CNT_CORR", "Protein_position", "Amino_acids", "Codons",
						"HGVSp", "Existing_variation",
						"Consequence", "IMPACT", "SIFT", "PolyPhen",
						"HGVSc", "CHROM", "POS", "MANE_SELECT","Gene", "EXON", "cDNA_position",	"CDS_position",	 
						"QUAL", "BASE_QUAL", "FILTER", "TYPE", "CTX",
						"SB_P", "AA_SKEW", "HP_FRAC",	"HP_LEN")
	} else {
		vcf_info_df = as.data.frame(vcf_geno[[1]])
		for (i in seq_along(vcf_geno)) {
			vcf_info_df[[i]] <- vcf_geno[[i]]
		}
		colnames(vcf_info_df) <- names(vcf_geno)
		keep_columns <- c(names(vcf_geno), "SYMBOL", "Protein_position", "Amino_acids", "Codons", "HGVSp", "Existing_variation", "Consequence", "IMPACT", "SIFT", "PolyPhen", "HGVSc", "MANE_SELECT", "Gene", "EXON", "cDNA_position", "CDS_position")
		cli_alert_info("Using data from FILTER and SAMPLE instead of INFO")
	}
	# selecting fixed columns
	CSQ_colnames <- c("Allele","Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
						"Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
						"cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", 
						"MANE_SELECT", "SIFT", "PolyPhen")

	CSQ <- ensemblVEP::parseCSQToGRanges(vcf, VCFRowID=rownames(vcf))

	ref <- str_split(rownames(meta(vcf_header)['contig'][[1]])[1], "\\.")[[1]][1]
	print(ref)
	CSQ_filt <- CSQ[!is.na(CSQ$Gene) & CSQ$Feature == ref]  # keeping only gene annotation and entries of the reference of intererst;

	vcf_CSQ_df_raw <- BiocGenerics::as.data.frame(CSQ_filt@elementMetadata) 
	rownames(vcf_CSQ_df_raw) <- names(CSQ_filt)[vcf_CSQ_df_raw$VCFRowID] # using CSQ_filt_df$VCFRowID assigned before to extract unique variant names

	cli_alert_info("Missing selected columns in CSQ: {sum(!(CSQ_colnames %in% colnames(vcf_CSQ_df_raw)))} columns.")

	vcf_CSQ_df <- vcf_CSQ_df_raw %>%
		select(all_of(CSQ_colnames))

	# Alternative 1:
	# dplyr::mutate(CSQ = as.character(CSQ)) %>%
	# tidyr::separate(col=CSQ, into=vcf_header_info_CSQ_colnames, sep="\\|")

	#Alternative 2:
	# [ ] - see how processing is done inside the ensemblVEP::parseCSQToGRanges and use directly this without calling ensemblVEP package!
	#ensemblVEP::parseCSQToGRanges
	# https://github.com/Bioconductor/ensemblVEP/blob/master/R/methods-parseCSQToGRanges.R

	# rename some of the columns to match Niko's names and re-order; remove empty annotation columns!
	# save as excel and tsv or csv
	# make into production and process flongle20
	# fix duplicated names
	cli_progress_step("Generating combined table.") # Removing empty columns.
	# to make sure if rownames == ID in vcf_ranges_df
	cli_alert_info("Check rownames equal ID: {all(rownames(vcf_ranges_df) == vcf_ranges_df$ID)} ")
	rownames(vcf_ranges_df) <- vcf_ranges_df$ID  # assigning rownames as these were assigned to ID at the beginning
	cli_alert_info("Check variant names before combining - identical: {all(rownames(vcf_ranges_df) == rownames(vcf_info_df), rownames(vcf_ranges_df) == rownames(vcf_CSQ_df))} ")

	# to make sure rownames are matched
	vcf_info_df$ID <- rownames(vcf_info_df)
	vcf_CSQ_df$ID <- rownames(vcf_CSQ_df)

	print(vcf_info_df)
	print(vcf_CSQ_df)

	vcf_complete_df_raw <- vcf_ranges_df %>%
		left_join(., vcf_info_df, by = "ID") %>%
		left_join(., vcf_CSQ_df, by = "ID")

	cli_alert_info("duplicated names in the final vcf: {sum(duplicated(names(vcf_complete_df_raw)))}")

	

	vcf_complete_df <- vcf_complete_df_raw %>%
	select(all_of(keep_columns))

	# saving as tsv and excel ----
	# [ ] - add check if the file was correctly processed
	output_tsv <- paste0(output_basename, ".tsv")
	output_excel <- paste0(output_basename, ".xlsx")
	cli_progress_step("Saving table as tsv and excel file.")
	cli_alert_info("tsv file: {output_tsv}")
	cli_alert_info("excel file: {output_excel}")

	vcf_complete_df <- data.frame(lapply(vcf_complete_df, as.character), stringsAsFactors=FALSE)
	write_tsv(vcf_complete_df, file = output_tsv)
	write.xlsx(vcf_complete_df, file = output_excel, asTable = TRUE, overwrite = TRUE)

	#cli_status_clear(id = sb)
	cli_alert_success("Vcf successfully processed.")
	cli_progress_done()
}