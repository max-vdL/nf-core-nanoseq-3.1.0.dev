#!/usr/bin/env Rscript

# script to convert from reference to genomic coordinates 
# using ensembl.db package

# TO-DO
# ! normalize vcf! (atomize?) - check if this is done in Niko's script!?
# check how the vcf corresponds to the tsv produced by Niko! - number of SNVs
# - [ ] add ensemblDB only for individual chromosomes by matchin the ENST-ID (e.g. chr9 for ABL1)
# - [ ] add asserts and test!
# - [x] speedup vcf loading part - load only coordinates for transformation and nothing else?
# - [ ] ? paralellize transcriptToGenome
# - [ ] speed up by rewriting in python - https://bioinformatics.stackexchange.com/questions/16353/ensembldb-equivalent-in-python
# - [ ] https://blog.liang2.tw/posts/2017/11/use-ensdb-database-in-python/
# - [ ] library+explicit(::) vs import??? 

# parallelization:
# - purrr, furrr, lapply and parallel versions

# - [ ] check warning from writing vcf
# Warning message:
#   In .Call(.make_vcf_geno, filename, fixed, names(geno), as.list(geno),  :
#              converting NULL pointer to R NULL

# importing functions
import::from(.from = optparse, make_option, OptionParser, parse_args, print_help)
import::from(.from = cli, cli_rule, cli_progress_bar, cli_alert_info, cli_alert_warning, cli_progress_step, cat_line, cli_alert_success, cli_progress_done)


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="vcf file in transcript coordinates - *.filt.vcf", metavar="character"),
  make_option(c("-c", "--cache"), type="character", default=".cache/AnnotationHub", 
              help="specify path to AnnotationHub cache", metavar="character"),
  make_option(c("--transcript"), type="character", default=NULL,
              help="specify the Ensembl Transcript ID that will be used to map coordinates")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# checking input files
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Input vcf file needs to be supplied!", call.=FALSE)
}
if (is.null(opt$transcript)){
  print_help(opt_parser)
  stop("Ensembl Transcript ID needs to be supplied!", call.=FALSE)
}

# check if vcf file exists
input_vcf <- file.path(opt$file)
if (!file.exists(input_vcf)){
  stop("Vcf file does not exists!", call.=FALSE)
}

# importing additional functions ----
# not needed to import if the above checks fail
import::from(.from = AnnotationHub, setAnnotationHubOption, AnnotationHub, snapshotDate, query)
import::from(.from = ensembldb, EnsDb, filter, transcriptToGenome)
import::from(.from = VariantAnnotation, readVcf, isSNV, isIndel, writeVcf)
import::from(.from = SummarizedExperiment, rowRanges)
import::from(.from = IRanges, ranges)

cli_rule(left = "transcript2genome_coord")
cli_progress_bar("Converting coordinates", type="tasks")

# specify fixed parameters
cache_dir = opt$cache
cli_alert_info("cache_dir: {cache_dir}")
AnnotationHub::setAnnotationHubOption("CACHE", cache_dir)

# save database only for chromosome9!
# save in the cache directory otherwise generate

ensembl_version=103
cache_db_file <- file.path(paste0(cache_dir, "/EnsDb_Hsapiens_v",ensembl_version,"_chr9.sqlite"))

if (!file.exists(cache_db_file)) {
  cli_alert_warning("AnnotationHub cache not detected.")
  cli_progress_step("Generating cache now.")
  cat_line()  # added to to print snapshotDate(): 2020-10-27 on the new line
  
  ah <- AnnotationHub::AnnotationHub()
  AnnotationHub::snapshotDate(ah) <- "2020-10-27"
  ahDb <- AnnotationHub::query(ah, pattern = c("Homo Sapiens", "EnsDb", ensembl_version))
  ahEdb <- ahDb[["AH89426"]]
  
  # only chr9 where ABL1 is
  edb_chr9 <- ensembldb::filter(ahEdb, filter = ~ seq_name == "9")
  
  # change below to ahEdb if all chromosomes
  old_cache_db_file <- paste0(edb_chr9@ensdb@dbname)
  
  file.copy(from = old_cache_db_file, to = cache_db_file)

  ahEdb <- ensembldb::EnsDb(cache_db_file)
} else {
  cli_progress_step("Loading cache now.")
  ahEdb <- ensembldb::EnsDb(cache_db_file)
}

# create output vcf name from input vcf name
output_vcf <- gsub(pattern = ".filt.vcf", replacement = ".filt.gmap.vcf", x=input_vcf)

# check if vcf is empty
cli_progress_step("Loading vcf file.")
vcf_read <- readLines(input_vcf)
vcf_read <- vcf_read[!grepl("#", vcf_read)]
if(length(vcf_read)==0){
  cli_alert_success("Empty vcf file cannot be mapped. Saving empty vcf file.")
  file.copy(input_vcf, output_vcf)
  cli_progress_done()
} else {
  # map vcf
  vcf <- VariantAnnotation::readVcf(input_vcf, genome = "GRCh38")
  cli_alert_info("Variants: {nrow(vcf)} SNVs: {sum(VariantAnnotation::isSNV(vcf))} INDELs: {sum(VariantAnnotation::isIndel(vcf))}")

  vcf_granges <- SummarizedExperiment::rowRanges(vcf)  # vcf@rowRanges
  vcf_iranges <- IRanges::ranges(vcf_granges) # vcf_granges@ranges

  # fix iranges names
  enst <- gsub("\\..*","",opt$transcript)
  iranges_names <- names(vcf_iranges)
  iranges_names_correct <- gsub(pattern = paste0("(", enst, ")(.+)"), replacement = "\\1", x = iranges_names)
  names(vcf_iranges) <- iranges_names_correct

  # production ----
  # [ ] - add tests to check outputs!
  # [ ] - add progress message
  cli_progress_step(msg = "Mapping coordinates.")
  rng_gnm <- ensembldb::transcriptToGenome(vcf_iranges, ahEdb) # - [ ] does my transcript definition contain 5' UTR (transcriptToGenome assumes including the 5' UTR)
  #rng_prt <- transcriptToProtein(vcf_iranges, ahEdb)

  cli_alert_info("Number of entries before and after conversion equals: {all(length(vcf_iranges) == length(rng_gnm))}")

  rng_granges_gnm <- unlist(rng_gnm)
  names(rng_granges_gnm) <- iranges_names

  # replace ranges
  vcf@rowRanges <- rng_granges_gnm

  cli_progress_step("Saving converted vcf file.")
  cli_alert_info("New vcf file: {output_vcf}")
  VariantAnnotation::writeVcf(obj = vcf, 
          file = output_vcf, 
          index=FALSE)

  cli_alert_success("Vcf file successfully processed.")
  cli_progress_done()
}