suppressWarnings(suppressPackageStartupMessages({
    library(glue)
    library(gwasvcf)
    library(VariantAnnotation)
    library(dplyr)
    library(magrittr)
    library(TwoSampleMR)
    library(MRInstruments)
    library(gwasglue)
    library(data.table)
    library(ieugwasr)
    library(plyr)
    library(mrpipeline)
    library(tidyr)
    library(xlsx)
    library(writexl)
    library(argparse)
    library(unixtools)
    #library("doParallel")
    library(parallel)
    library(MendelianRandomization)
}))

tag_r2=0.8
n_vcf_parallel=10
N=10 #at a time how many exposure should run

set_bcftools("/usr/local/bin/bcftools")
set_plink("/usr/local/bin/plink")
# Create argument parser
parser <- ArgumentParser()

# Add arguments
parser$add_argument("--exposure_file",
                    required = TRUE,
                    help = "Details of the exposure GWAS. The file should contain two columns: 'gwas_id' and 'file_path'. The 'file_path' column should contain the full path to the GWAS summary statistics in VCF format.; multi allele split to bi allele")

parser$add_argument("--outcome_file",
                    required = TRUE,
                    help = "Details of the outcome GWAS. The file should contain two columns: 'gwas_id' and 'file_path'. The 'file_path' column should contain the full path to the GWAS summary statistics in VCF format.; multi allele split to bi allele")

parser$add_argument("--output_folder",
                    required = TRUE,
                    help = "Name of the folder where the MR analysis results should be saved.")

parser$add_argument("--ld_file",
                    required = TRUE,
                    help = "File containing the linkage disequilibrium data.,file prefix with full path")

parser$add_argument("--variants_include",
                    default = "NA",
                    help = "A file containing variant IDs that should be include ; all these variants should be present in the ld_file")

parser$add_argument("--infoscore",
                    help = "Info score cutoff to remove low-quality variants.",
                    default = 0,
                    type = "double")

parser$add_argument("--maf",
                    help = "MAF cutoff to remove variants with low minor allele frequency.",
                    type = "double",
                    default = 1)

parser$add_argument("--pvalue",
                    help = "Variants with higher than this p value will not be included in the exposure vcf file",
                    type = "double",
                    default = 1)

parser$add_argument("--indels",
                    help = "Specify whether insertion deletion variations should be considered ('Yes' or 'No').",
                    choices = c("Yes", "No"),
                    default = "No")

parser$add_argument("--variant_id",
                    help = "Specify whether variant id is rsid or unique_id.",
                    choices = c("rsid", "unique_id"),
                    default = "unique_id")

parser$add_argument("--palindromic",
                    help = "Specify whether palindromic variations should be considered ('Yes' or 'No').",
                    choices = c("Yes", "No"),
                    default = "No")

parser$add_argument("--dbsnpfile",
                    help = "file with rsid (full path should be provided); it should contain five columns , namely 'seqnames','start','ID','REF','ALT' and the columns are in the same order as mentioned",
                    default = "NA")

parser$add_argument("--mhc_remove",
                    help = "Specify whether to remove MHC region or not",
                    choices = c("Yes", "No"),
                    default = "No")

parser$add_argument("--mhc_region",
                    help = "Mhc region to be removed",
                    choices = c("6:25000000-35000000", "6:28477797-33448354", "6:28510120-33480577"),
                    default = "6:25000000-35000000")

parser$add_argument("--remove_chrX",
                    help = "Specify whether the variants in chromosome X should be removed or not",
                    choices = c("Yes", "No"),
                    default = "No")

# exposure_file="/mnt/disks/sdd/ukb_b_sumstat_vcf/ieu-b_mr/exposure_sumstat_vcf_files.csv"
# outcome_file="/mnt/disks/sdd/ukb_b_sumstat_vcf/ieu-b_mr/outcome_sumstat.csv"
# output_folder="/mnt/disks/sdd/ukb_b_sumstat_vcf/ieu-b_mr/ieu-b_mr2/"
# ld_file="/mnt/disks/sdd/test_pipelines/onekg_ref_files/onekg_vcf_files/plink_files/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_multiallele_uniqid_SNP"
# variants_include="/mnt/disks/sdd/test_pipelines/onekg_ref_files/onekg_vcf_files/plink_files/EUR_Variants.txt"
# infoscore_cutoff=0.6
# maf_cutoff=0.01
# pvalue_cutoff=0.00000005
# indels="Yes"
# variant_id_type="unique_id"
# palindromic="Yes"
# dbsnpfile="NA"
# mhc_remove="Yes"
# mhc_region="6:25000000-35000000"
# remove_chrX="Yes"



# Parse arguments
args <- parser$parse_args()

# Retrieve argument values
exposure_file <- args$exposure_file
outcome_file <- args$outcome_file
output_folder <- args$output_folder
ld_file <- args$ld_file
variants_include <- args$variants_include
infoscore_cutoff <- args$infoscore
maf_cutoff <- args$maf
pvalue_cutoff <- args$pvalue
indels <- args$indels
variant_id_type <- args$variant_id
palindromic <- args$palindromic
dbsnpfile <- args$dbsnpfile
mhc_remove <- args$mhc_remove
mhc_region <- args$mhc_region
remove_chrX <- args$remove_chrX

exposure_file_df=fread(exposure_file)
outcome_file_df=fread(outcome_file)
lower_AF=maf_cutoff
upper_AF=1-maf_cutoff
tmp_folder=glue("{output_folder}/tmp_folder")
system(glue("mkdir -p {tmp_folder}"))
set.tempdir(tmp_folder)


print("Started the Analysis")

# Define the vcf_filter function
vcf_filter <- function(i, outcome_file_df, lower_AF, upper_AF, infoscore_cutoff, filtered_vcf_folder, variant_id_type = "unique_id", p_value = 1,variants_include="NA") {
    filt_vcf_folder <- glue("{output_folder}/{filtered_vcf_folder}")
    dir.create(filt_vcf_folder, showWarnings = FALSE, recursive = TRUE)
    
    sample_id <- as.character(outcome_file_df[i, 'gwas_id'])
    vcf_file <- as.character(outcome_file_df[i, 'file_path'])
    
    LP_value <- -log10(p_value)
    
    # Get the VCF format
    vcformat <- strsplit(system(glue("zgrep -v '#' {vcf_file} | head -n1 | awk '{{print $9}}'"), intern = TRUE), ":")[[1]]
    
    indels <- "No"  # Assuming indels is fixed; adjust as needed or pass it as a parameter
    conditions <- c()
    
    # Build conditions based on available format fields
    if ("AF" %in% vcformat) {
        conditions <- c(conditions, glue("FORMAT/AF>{lower_AF} & FORMAT/AF<{upper_AF}"))
    }
    if ("SI" %in% vcformat) {
        conditions <- c(conditions, glue("FORMAT/SI>{infoscore_cutoff}"))
    }
    
    if ("LP" %in% vcformat) {
        conditions <- c(conditions, glue("FORMAT/LP>={p_value}"))
    }
    

    if (variants_include!="NA") {
        conditions <- c(conditions, glue(" 'ID=@{variants_include}' "))
    }

    id="'%CHROM\\_%POS\\_%REF\\_%ALT'"
    # If there are conditions to filter by
    if (length(conditions) > 0) {
        filter_condition <- paste(conditions, collapse = " & ")
        vcf_file_out <- glue("{filt_vcf_folder}/{sample_id}_filt.vcf.gz")
        
        # Build the command based on the value of indels and variant_id_type
        if (tolower(indels) == "no") {
            if (variant_id_type == "unique_id") {
                command <- glue(" bcftools annotate --set-id {id}  {vcf_file}  | bcftools view --types snps --include '{filter_condition}' | bgzip -c > {vcf_file_out}")
            } else {
                command <- glue("bcftools view --types snps --include '{filter_condition}' {vcf_file} | bgzip -c > {vcf_file_out}")
            }
        } else {
            if (variant_id_type == "unique_id") {
                command <- glue(" bcftools annotate --set-id {id} {vcf_file}  | bcftools view --include '{filter_condition}' {vcf_file} | bgzip -c > {vcf_file_out}")
            } else {
                command <- glue("bcftools view --include '{filter_condition}' {vcf_file} | bgzip -c > {vcf_file_out}")
            }
        }
        
        # Execute the command
        system(command) 
        # Execute tabix command with error handling
        tabix_command <- glue("tabix -f -p vcf {vcf_file_out}")
        tryCatch({
            system(tabix_command)
        }, error = function(e) {
            message(glue("Error executing tabix command for {vcf_file_out}: {e$message}"))
        })
    }
}


##GWAS VCF TO OUTCOME DATA ; during proxy search this was not working perfectly "snp_col" value required "SNP"
gwasvcf_to_TwoSampleMR2<-function (vcf, type = "exposure") 
        {a <- vcf %>% gwasvcf::vcf_to_granges()
            S4Vectors::mcols(a)[["SNP"]] <- names(a)
            a <- dplyr::as_tibble(a)
            if (!"ES" %in% names(a)) 
                a[["ES"]] <- NA
            if (!"SE" %in% names(a)) 
                a[["SE"]] <- NA
            if (!"LP" %in% names(a)) 
                a[["LP"]] <- NA
            if (!"SS" %in% names(a)) 
                a[["SS"]] <- NA
            if (!"NC" %in% names(a)) 
                a[["NC"]] <- NA
            if (!"id" %in% names(a)) 
                a[["id"]] <- NA
            a[["LP"]] <- 10^-a[["LP"]]
            a[["NCONT"]] <- a[["SS"]] - a[["NC"]]
            TwoSampleMR::format_data(a, type = type, snp_col = "ID", 
                effect_allele_col = "ALT", other_allele_col = "REF", 
                eaf_col = "AF", chr_col = "seqnames", pos_col = "start", 
                beta_col = "ES", se_col = "SE", pval_col = "LP", samplesize_col = "SS", 
                ncase_col = "NC", ncontrol_col = "NCONT", phenotype_col = "id")
        }



# LD pruning function
ld_pruning <- function(exposure_id, df2, ld_file, local_qc_passed_exposure_variants_failed_during_ld_prune_df, local_qc_passed_exposure_variants_after_ld_prune_df) {
    df3 <- data.frame()

    # Perform LD clumping and handle potential errors
    tryCatch({
        df3 <- suppressMessages(suppressWarnings({
            ld_clump(dplyr::tibble(rsid = df2$ID, pval = df2$P, id = df2$id),
                    plink_bin = "/usr/local/bin/plink",
                    bfile = ld_file)
        }))
    }, error = function(e) {
        # Log the error if needed
        message(paste("LD clumping failed for", exposure_id, ":", e$message))
    })

    # Check if LD clumping resulted in an empty data frame
    if (nrow(df3) == 0) {
        selected_df <- df2
        if (nrow(selected_df) > 0) {
            tmp_selected_df <- data.table::copy(selected_df)
            tmp_selected_df$exposure_id <- exposure_id
            local_qc_passed_exposure_variants_failed_during_ld_prune_df <- dplyr::bind_rows(local_qc_passed_exposure_variants_failed_during_ld_prune_df, tmp_selected_df)
        }
    } else {
        # Ensure df3 has the correct columns and merge with df2
        df3 <- as.data.frame(df3)
        colnames(df3) <- c("ID", 'pval',"id")
        selected_df <- merge(df2, df3, by = c("ID", "id"), all = FALSE)
        if (nrow(selected_df) > 0) {
            tmp_selected_df <- data.table::copy(selected_df)
            tmp_selected_df$exposure_id <- exposure_id
            local_qc_passed_exposure_variants_after_ld_prune_df <- dplyr::bind_rows(local_qc_passed_exposure_variants_after_ld_prune_df, tmp_selected_df)
        }
    }
    
    result <- list(
        selected_df = selected_df,
        local_qc_passed_exposure_variants_failed_during_ld_prune_df = local_qc_passed_exposure_variants_failed_during_ld_prune_df,
        local_qc_passed_exposure_variants_after_ld_prune_df = local_qc_passed_exposure_variants_after_ld_prune_df
    )
    return(result)
}


# Function to read, process exposure VCF, and perform LD pruning
read_exposure_vcf_and_process <- function(exposure_id,dbsnp_df, ld_file, 
                                          infoscore_cutoff, maf_cutoff, pvalue_cutoff, variant_id_type,
                                          mhc_remove, mhc_region, remove_chrX,
                                          local_qc_passed_exposure_variants_before_ld_prune_df,
                                          local_qc_passed_exposure_variants_after_ld_prune_df,
                                          local_qc_passed_exposure_variants_failed_during_ld_prune_df,
                                          local_exposure_with_no_significant_variants_before_ld_prune,
                                          local_exposure_with_no_significant_variants_after_ld_prune) {
    
    df <- data.frame()
    exposure_df <- data.frame()
    
    # Read VCF file
    tryCatch({
        vcf <- readVcf(glue("{output_folder}filtered_exposure_vcf/{exposure_id}_filt.vcf.gz"))
        df <- vcf_to_granges(vcf) %>% as_tibble()
    }, error = function(e) {
        message(glue("Error reading VCF file for sample '{exposure_id}': {conditionMessage(e)}"))
    })
    
    if (nrow(df) > 0) {
        df<-as.data.frame(df)
        # Convert specified columns to numeric
        cols_to_convert <- c("SI", "AF", "LP", "start", "ES", "SE")
        df[cols_to_convert] <- lapply(df[cols_to_convert], as.numeric)
        
        # Merge with dbSNP data if provided
        if (dbsnp_df!="NA") {
            df <- merge(df, dbsnp_df, by = c('seqnames', 'start', 'REF', 'ALT'))
        } else {
            message("dbSNP data not provided")
        }
        
        # Add unique ID to variant if required
        if (variant_id_type == 'unique_id') {
            df$ID <- paste(df$seqnames, df$start, df$REF, df$ALT, sep = "_")
        }
        
        # Filter variants based on Infoscore
        if ("SI" %in% colnames(df)) {
            df <- df[df$SI >= infoscore_cutoff | is.na(df$SI), ]
            message(glue("Filtering variants based on infoscore for sample '{exposure_id}': infoscore cutoff was {infoscore_cutoff}"))
        } else {
            message(glue("Imputation score column not present or contains only missing values for sample '{exposure_id}', so filtering based on Imputation score not performed"))
        }
        
        # Filter variants based on MAF
        if ("AF" %in% colnames(df)) {
            df$MAF <- ifelse(df$AF > 0.5, 1 - df$AF, df$AF)
            df <- df[df$MAF >= maf_cutoff | is.na(df$AF), ]
            message(glue("Filtering variants based on MAF for sample '{exposure_id}': MAF cutoff was {maf_cutoff}"))
        } else {
            message(glue("MAF column not present or contains only missing values for sample '{exposure_id}', so filtering based on MAF not performed"))
        }
        
        # Remove Palindromic variants if MAF > 0.42
        if (palindromic == 'Yes' && "AF" %in% colnames(df)) {
            condition <- (df$MAF <= 0.42 | is.na(df$MAF)) | !(((df$REF == "C" & df$ALT == "G") | (df$REF == "G" & df$ALT == "C")) | 
                                                                   ((df$REF == "A" & df$ALT == "T") | (df$REF == "T" & df$ALT == "A")))
            df <- df[condition, ]
            message(glue("Removing palindromic variants if MAF > 0.42 for sample '{exposure_id}'"))
        } else {
            message(glue("Palindromic variants not removed from sample '{exposure_id}'"))
        }
        
        # Filter based on P value
        if ("LP" %in% colnames(df)) {
            df$P <- 10 ^ -df$LP
            df <- df[df$P <= pvalue_cutoff | is.na(df$LP), ]
            message(glue("Filtering variants based on pvalue cutoff for sample '{exposure_id}': pvalue cutoff was {pvalue_cutoff}"))
        } else {
            message(glue("LP column not present or contains only missing values for sample '{exposure_id}', so filtering based on pvalue not performed"))
        }
        
        # Filter based on MHC regions
        if (mhc_remove == "Yes" && !is.null(mhc_region)) {
            mhc_region <- strsplit(mhc_region, "[:-]")[[1]]
            mhc_chr <- as.character(mhc_region[1])
            mhc_start <- as.integer(mhc_region[2])
            mhc_end <- as.integer(mhc_region[3])
            df <- df[!(as.character(df$seqnames) == mhc_chr & df$start > mhc_start & df$start < mhc_end), ]
            message(glue("Removing variants present in the MHC region for sample '{exposure_id}': MHC Region was {mhc_chr}:{mhc_start}-{mhc_end}"))
        } else {
            message(glue("Not removing variants from the MHC regions for sample '{exposure_id}'"))
        }
        
        # Remove variants from X chromosome
        if (remove_chrX == 'Yes') {
            df <- df[as.character(df$seqnames) != "X", ]
            message(glue("Removing variants present in X chromosome for sample '{exposure_id}'"))
        } else {
            message(glue("Variants present in X chromosome not removed from sample '{exposure_id}'"))
        }
        
        # Remove variants from Y chromosome and MT
        df <- df[!(as.character(df$seqnames) %in% c("Y", "MT")), ]
        message(glue("Removing variants present in Y and MT chromosome for sample '{exposure_id}'"))
        
        # Perform LD pruning
        if (nrow(df) > 0) {
            tmp_df <- df
            tmp_df$exposure_id <- exposure_id
            local_qc_passed_exposure_variants_before_ld_prune_df <- dplyr::bind_rows(local_qc_passed_exposure_variants_before_ld_prune_df, tmp_df)
            message(glue("Total number of QC passed variants before LD pruning for the sample '{exposure_id}' is {nrow(df)}"))
            
            if (nrow(df) > 1) {
                tryCatch({
                    result <- ld_pruning(exposure_id, df, ld_file, local_qc_passed_exposure_variants_failed_during_ld_prune_df, local_qc_passed_exposure_variants_after_ld_prune_df)
                    exposure_df<-result$selected_df
                    local_qc_passed_exposure_variants_failed_during_ld_prune_df<-result$local_qc_passed_exposure_variants_failed_during_ld_prune_df
                    local_qc_passed_exposure_variants_after_ld_prune_df<-result$local_qc_passed_exposure_variants_after_ld_prune_df

                }, error = function(e) {
                    message(glue("LD pruning failed for '{exposure_id}', using all the exposure variants: {conditionMessage(e)}"))
                    exposure_df <- df
                })
                
            } else {
                exposure_df <- df
            }
        } else {
            message(glue("No QC passed variants for sample '{exposure_id}'"))
            local_exposure_with_no_significant_variants_before_ld_prune <- c(local_exposure_with_no_significant_variants_before_ld_prune, exposure_id)
        }
    } else {
        message(glue("No QC passed variants for sample '{exposure_id}'"))
        local_exposure_with_no_significant_variants_after_ld_prune <- c(local_exposure_with_no_significant_variants_after_ld_prune, exposure_id)
    }
    message(glue("Total number of QC passed variants after LD pruning for the sample '{exposure_id}' is {nrow(exposure_df)}"))
    results <- list(
        exposure_df = exposure_df,
        local_qc_passed_exposure_variants_before_ld_prune_df = local_qc_passed_exposure_variants_before_ld_prune_df,
        local_qc_passed_exposure_variants_after_ld_prune_df = local_qc_passed_exposure_variants_after_ld_prune_df,
        local_qc_passed_exposure_variants_failed_during_ld_prune_df = local_qc_passed_exposure_variants_failed_during_ld_prune_df,
        local_exposure_with_no_significant_variants_before_ld_prune = local_exposure_with_no_significant_variants_before_ld_prune,
        local_exposure_with_no_significant_variants_after_ld_prune = local_exposure_with_no_significant_variants_after_ld_prune
        )
    return(results)
}


###################
format_exposure <- function(exposure_id,selecteed_df,pvalue_cutoff) {
                    exp_df <- format_data(
                    selecteed_df,
                    type = "exposure",snps = NULL,header = TRUE,
                    phenotype_col = "exposure",snp_col = "ID",
                    beta_col = "ES",se_col = "SE",
                    eaf_col = "AF",
                    effect_allele_col = "ALT",other_allele_col = "REF",
                    pval_col = "P",
                    samplesize_col = "SS",
                    gene_col = "id",min_pval = pvalue_cutoff,
                    chr_col = "seqnames",pos_col = "start",
                    log_pval = FALSE)
            return(exp_df) }


#format Outcome data
format_outcome <- function(outcome_id, exposure_id, outcomevcf, formatted_selected_exposure_df, ld_file,tag_r2) {
    # First try-catch for query_gwas
    vcf <- tryCatch({
        query_gwas(outcomevcf, rsid = formatted_selected_exposure_df$SNP, proxies = "yes", bfile = ld_file, tag_r2 =tag_r2)
    }, error = function(e) {
        message(glue("Failed to query GWAS for outcome_id: {outcome_id} and exposure_id: {exposure_id}. Error: {e$message}"))
        return(NULL)
    })

    if ( nrow(as.data.frame(vcf_to_granges(vcf) %>% as_tibble()))<1 ) {
        local_outcome_with_no_variant_df <- rbind(local_outcome_with_no_variant_df, data.frame(outcome_id = outcome_id, exposure_id = exposure_id, stringsAsFactors = FALSE))
    }

    # Second try-catch for gwasvcf_to_TwoSampleMR
    outcome_df <- tryCatch({
        gwasvcf_to_TwoSampleMR(vcf, type = "outcome")
    }, error = function(e) {
        message(glue("First method for conversion of query_gwas output from outcome GWAS: {outcome_id} to dataframe failed for outcome_id: {outcome_id} and exposure_id: {exposure_id}."))
        return(NULL)
    })

    if (is.null(outcome_df)) {
        # Third try-catch for gwasvcf_to_TwoSampleMR2
        outcome_df <- tryCatch({
            gwasvcf_to_TwoSampleMR2(vcf, type = "outcome")
        }, error = function(e) {
            message(glue("Second method for conversion of query_gwas output from outcome GWAS: {outcome_id} to dataframe failed for outcome_id: {outcome_id} and exposure_id: {exposure_id}."))
            local_outcome_with_no_variant_df <- rbind(local_outcome_with_no_variant_df, data.frame(outcome_id = outcome_id, exposure_id = exposure_id, stringsAsFactors = FALSE))
            return(data.frame())
        })
    }

    return(outcome_df)
}


###This is to harmonise exposure and outcome data from MR Analysis input
Harmonise_exposure_Outcome_data <- function(outcome_id, exposure_id, exposure_df, outcomevcf, 
                                                        ld_file, pvalue_cutoff,tag_r2,
                                                        local_outcome_with_no_variant_df,
                                                        local_harmonised_df){

    dat <- data.frame()
    formatted_selected_exposure_df <- data.frame()

    # Remove column 'V1' if it exists
    if ("V1" %in% colnames(exposure_df)) {
        exposure_df <- exposure_df[, !names(exposure_df) %in% "V1", with = FALSE]
    }

    # Add an 'exposure' column
    exposure_df$exposure <- exposure_df$id

    # Select the unique exposure
    Exposure <- sort(unique(exposure_df$exposure))


    # Check if the filtered exposure data frame is not empty
    if (nrow(exposure_df) > 0) {
        formatted_selected_exposure_df <- format_exposure(exposure_id,exposure_df, pvalue_cutoff)

            if (variant_id_type == 'unique_id') {
                formatted_selected_exposure_df$SNP<-toupper(formatted_selected_exposure_df$SNP) 
            }

        if (nrow(formatted_selected_exposure_df) > 0) {
            formatted_selected_exposure_df <- calc_f_stat(formatted_selected_exposure_df, f_cutoff = 0)
            formatted_selected_exposure_df$exposure <- Exposure

            # Try to format the outcome data
            formatted_outcome_df <- tryCatch({
            format_outcome(outcome_id, exposure_id, outcomevcf, formatted_selected_exposure_df, ld_file,tag_r2)
            }, error = function(e) {
                message(glue("No SNPs present in the outcome VCF file for outcome_id: {outcome_id} and exposure_id: {exposure_id}. Error: {e$message}"))
                return(data.frame())
            })
        }
    }

        # Check if the formatted outcome data frame is not empty before harmonizing the data
    if (nrow(formatted_outcome_df) > 0) {
            if (variant_id_type == 'unique_id') {
                formatted_outcome_df$SNP<-toupper(formatted_outcome_df$SNP) 
            }

        dat <- tryCatch({
            dat <- harmonise_data(exposure_dat = formatted_selected_exposure_df, outcome_dat = formatted_outcome_df, action = 1)
        }, error = function(e) {
            message(glue("Error during harmonization for outcome_id: {outcome_id} and exposure_id: {exposure_id}. Error: {e$message}"))
            return(data.frame())
        })
    }

    if (nrow(dat) > 0) {
        local_harmonised_df <- rbind.fill(local_harmonised_df, dat)
    }

    if (nrow(formatted_outcome_df)<1) {
        local_outcome_with_no_variant_df <- rbind(local_outcome_with_no_variant_df, data.frame(outcome_id = outcome_id, exposure_id = exposure_id, stringsAsFactors = FALSE))
    }
    results <- list(
        dat = dat,
        local_outcome_with_no_variant_df = local_outcome_with_no_variant_df,
        local_harmonised_df = local_harmonised_df)

        return(results)
}


###MR Analysis
perform_two_sample_mr_analysis <- function(dat, outcome_id, exposure_id) {
    methodlist_1 <- c("mr_wald_ratio", "mr_two_sample_ml", "mr_egger_regression", "mr_egger_regression_bootstrap", "mr_simple_median", "mr_weighted_median", "mr_penalised_weighted_median", "mr_ivw", "mr_ivw_radial", "mr_ivw_mre", "mr_ivw_fe", "mr_simple_mode", "mr_weighted_mode", "mr_weighted_mode_nome", "mr_simple_mode_nome", "mr_sign", "mr_uwr")
    methodlist_2 <- c("mr_wald_ratio", "mr_two_sample_ml", "mr_egger_regression", "mr_egger_regression_bootstrap", "mr_simple_median", "mr_weighted_median", "mr_penalised_weighted_median", "mr_ivw", "mr_ivw_mre", "mr_ivw_fe", "mr_simple_mode", "mr_weighted_mode", "mr_weighted_mode_nome", "mr_simple_mode_nome", "mr_sign", "mr_uwr")
    methodlist_3 <- c("mr_wald_ratio", "mr_two_sample_ml", "mr_egger_regression", "mr_egger_regression_bootstrap", "mr_simple_median", "mr_weighted_median", "mr_penalised_weighted_median", "mr_ivw")
    heterogeneity_list_1 <- c("mr_two_sample_ml", "mr_egger_regression", "mr_ivw", "mr_ivw_radial", "mr_uwr")
    heterogeneity_list_2 <- c("mr_two_sample_ml", "mr_egger_regression", "mr_ivw", "mr_uwr")
    heterogeneity_list_3 <- c("mr_egger_regression", "mr_ivw", "mr_uwr")
    heterogeneity_list_4 <- c("mr_egger_regression")

    mr_res <- data.frame()
    mr_het <- data.frame()
    mr_pleo <- data.frame()
    mr_direction <- data.frame()

    tryCatch({
        mr_res <- mr(dat, method_list = methodlist_1)
        mr_res <- generate_odds_ratios(mr_res)
    }, error = function(e) {
        message(glue("TwosampleMR analysis with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed for methodlist_1 \n"))
    })

    if (nrow(mr_res) == 0) {
        tryCatch({
            mr_res <- mr(dat, method_list = methodlist_2)
            mr_res <- generate_odds_ratios(mr_res)
        }, error = function(e) {
            message(glue("TwosampleMR analysis with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed for methodlist_2 \n"))
        })
    }

    if (nrow(mr_res) == 0) {
        tryCatch({
            mr_res <- mr(dat, method_list = methodlist_3)
            mr_res <- generate_odds_ratios(mr_res)
        }, error = function(e) {
            message(glue("TwosampleMR analysis with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed for methodlist_3 \n"))
        })
    }

    tryCatch({
        mr_het <- mr_heterogeneity(dat, method_list = heterogeneity_list_1)
    }, error = function(e) {
        message(glue("Heterogeneity analysis with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed for heterogeneity_list_1 \n"))
    })

    if (nrow(mr_het) == 0) {
        tryCatch({
            mr_het <- mr_heterogeneity(dat, method_list = heterogeneity_list_2)
        }, error = function(e) {
            message(glue("Heterogeneity analysis with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed for heterogeneity_list_2 \n"))
        })
    }

    if (nrow(mr_het) == 0) {
        tryCatch({
            mr_het <- mr_heterogeneity(dat, method_list = heterogeneity_list_3)
        }, error = function(e) {
            message(glue("Heterogeneity analysis with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed for heterogeneity_list_3 \n"))
        })
    }
    if (nrow(mr_het) == 0) {
        tryCatch({
            mr_het <- mr_heterogeneity(dat, method_list = heterogeneity_list_4)
        }, error = function(e) {
            message(glue("Heterogeneity analysis with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed for heterogeneity_list_4 \n"))
        })
    }

    tryCatch({
        mr_pleo <- mr_pleiotropy_test(dat)
    }, error = function(e) {
        message(glue("Pleiotropy test with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed"))
    })

    tryCatch({
        mr_direction <- directionality_test(dat)
    }, error = function(e) {
        message(glue("Directionality test with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed"))
    })

    return(list(mr_res = mr_res, mr_het = mr_het, mr_pleo = mr_pleo, mr_direction = mr_direction))
}




##Single Variant annalysis
perform_single_variant_mr_analysis <- function(dat,outcome_id,exposure_id) {
    mr_single_wald_ratio <- NULL
    mr_single_meta_fixed <- NULL

    tryCatch({
        mr_single_wald_ratio <- mr_singlesnp(dat, parameters = default_parameters(), single_method = "mr_wald_ratio")
    }, error = function(e) {
    message(glue("mr_single_wald_ratio test with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed"))
    })

    tryCatch({
        mr_single_meta_fixed <- mr_singlesnp(dat, parameters = default_parameters(), single_method = "mr_meta_fixed")
    }, error = function(e) {
    message(glue("mr_single_meta_fixed test with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed"))
    })

    return(list(mr_single_wald_ratio = mr_single_wald_ratio, mr_single_meta_fixed = mr_single_meta_fixed))
}

################-----------------------------------------------------------------------------#####################################
#MendelianRandomization  https://github.com/cran/MendelianRandomization Version 0.8
perform_mendelian_randomization <- function(dat,outcome_id,exposure_id) {
    M_Randomization_df <- data.frame()
    # Initialize dat2
    dat2 <- dat_to_MRInput(dat)

    tryCatch({
        # Try running all methods
        MRAllObject_all <- MendelianRandomization::mr_allmethods(dat2[[1]], method = "all")
        # If successful, convert results to data frame
        M_Randomization_df <- as.data.frame(MRAllObject_all@Values)
        }, error = function(e) {
        message(glue("MendelianRandomization (https://github.com/cran/MendelianRandomization) test (with all method) with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed"))
        })

    # If all methods failed or returned no results, try IVW method
    if (nrow(M_Randomization_df) == 0) {
        tryCatch({
            MRAllObject_all <- MendelianRandomization::mr_allmethods(dat2[[1]], method = "ivw")
            M_Randomization_df <- as.data.frame(MRAllObject_all@Values)
        }, error = function(e) {
        message(glue("MendelianRandomization (https://github.com/cran/MendelianRandomization) test (with ivw) with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed"))
        })
    }

    # If there are results, add additional columns
    if (nrow(M_Randomization_df) > 0) {
        M_Randomization_df$outcome <- rep(dat$outcome[1], nrow(M_Randomization_df))
        M_Randomization_df$exposure <- rep(dat$exposure[1], nrow(M_Randomization_df))
        M_Randomization_df$SNPs <- paste(unique(dat$SNP), collapse = ",")
    }

    return(M_Randomization_df)
}

#MendelianRandomization  https://github.com/cran/MendelianRandomization IVW Delta Version 0.8
perform_mendelian_randomization_ivw_delta <- function(dat,outcome_id,exposure_id) {
    dat2 <- dat_to_MRInput(dat)
    mr_ivw_delta_df <- data.frame()

    tryCatch({
        t <- MendelianRandomization::mr_ivw(dat2[[1]], weights = "delta", distribution = "normal")
    }, error = function(e) {
        message(glue("MendelianRandomization (https://github.com/cran/MendelianRandomization) test (with ivw delta) with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed"))
        t <- data.frame() # Initialize t as an empty data frame in case of an error
    })

    if (typeof(t) == "S4") {
        mr_ivw_delta_df <- rbind.fill(mr_ivw_delta_df, data.frame(
            Model = t@Model,
            Outcome = t@Outcome,
            Penalized = t@Penalized,
            Estimate = t@Estimate,
            CILower = t@CILower,
            Alpha = t@Alpha,
            SNPs = t@SNPs,
            Heter.Stat = paste(unique(t@Heter.Stat), collapse = ","),
            Exposure = t@Exposure,
            Robust = t@Robust,
            Correlation = t@Correlation,
            StdError = t@StdError,
            CIUpper = t@CIUpper,
            Pvalue = t@Pvalue,
            RSE = t@RSE
        ))
    }

    return(mr_ivw_delta_df)
}

##perform_mr_presso
perform_mr_presso <- function(dat,outcome_id,exposure_id) {
    mr_presso_df <- data.frame()
    presso <- data.frame()

    tryCatch({
        mr_presso <- run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)
    }, error = function(e) {
        message(glue("mr_presso test with outcome_id: {outcome_id} and exposure_id: {exposure_id} failed"))
        mr_presso <- list() # Initialize mr_presso as an empty list in case of an error
    })

    tryCatch({
        presso <- mr_presso[[1]]$`Main MR results`
    }, error = function(e) {
        presso <- data.frame() # Initialize presso as an empty data frame in case of an error
    })

    if (nrow(presso) > 0) {
        presso$RSSobs <- mr_presso[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs
        presso$GlobalTest_Pvalue <- mr_presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
        presso$exposure <- attributes(mr_presso)$exposure
        presso$outcome <- attributes(mr_presso)$outcome
        mr_presso_df <<- rbind.fill(mr_presso_df, presso)
    }

    return(mr_presso_df)
}


#c("samplesize.outcome",'ncase.outcome','ncontrol.outcome','samplesize.exposure')
mendelian_analysis_pipeline <- function(exposure_part_df, part, outcome_file_df, output_folder, dbsnp_df, ld_file, 
                                              infoscore_cutoff, maf_cutoff, pvalue_cutoff, variant_id_type,mhc_remove, 
                                              mhc_region, remove_chrX,outcomevcf_path, tag_r2) {

    tryCatch({
        summary_folder=glue("{output_folder}/summary_exposure_outcome/")
        dir.create(glue("{summary_folder}"), recursive = TRUE)
    }, error = function(e) { cat("Error creating directory:", conditionMessage(e), "\n") })
                                    
    file_suffixes <- c("Harmonised_Exposure_Outcome", 'MendelianPipeline_Test', 'MendelianRandomization_AllTest',
                       'MendelianRandomization_IVW_Delta_Test', 'TwoSampleMR_Analysis_Directionality_Test',
                       'TwoSampleMR_Analysis_Heterogeneity_Test', 'TwoSampleMR_Analysis_Hpleiotropy_Test',
                       'TwoSampleMR_Analysis_Multiple_MR_Test', 'TwoSampleMR_Analysis_SingleVariantMetafixed_Test',
                       'TwoSampleMR_Analysis_SingleVariantWald_Test')

    for (exposure_id in exposure_part_df$gwas_id) {
        # Initialize local data frames
        local_qc_passed_exposure_variants_before_ld_prune_df <- data.frame()
        local_qc_passed_exposure_variants_after_ld_prune_df <- data.frame()
        local_qc_passed_exposure_variants_failed_during_ld_prune_df <- data.frame()
        local_exposure_with_no_significant_variants_before_ld_prune <- character()
        local_exposure_with_no_significant_variants_after_ld_prune <- character()

        cat("\n\n\n---------------------------------------------------\n")
        cat(glue("Started running {exposure_id}\n"))
        filt_exposure_vcf <- glue("{output_folder}filtered_exposure_vcf/{exposure_id}_filt.vcf.gz")
        variant_count <- as.integer(system(glue('zgrep -v "#" {filt_exposure_vcf} | wc -l'), intern = TRUE))

        if (variant_count > 0) {
            out_folder <- glue("{output_folder}{exposure_id}/")

            # Create the output directory and handle errors
            tryCatch({
                dir.create(out_folder, recursive = TRUE)
            }, error = function(e) { cat("Error creating directory:", conditionMessage(e), "\n") })

            results <- read_exposure_vcf_and_process(exposure_id, dbsnp_df, ld_file, 
                                                     infoscore_cutoff, maf_cutoff, pvalue_cutoff, variant_id_type,
                                                     mhc_remove, mhc_region, remove_chrX,
                                                     local_qc_passed_exposure_variants_before_ld_prune_df,
                                                     local_qc_passed_exposure_variants_after_ld_prune_df,
                                                     local_qc_passed_exposure_variants_failed_during_ld_prune_df,
                                                     local_exposure_with_no_significant_variants_before_ld_prune,
                                                     local_exposure_with_no_significant_variants_after_ld_prune) 

            exposure_df <- results$exposure_df
            if (nrow(results$local_qc_passed_exposure_variants_before_ld_prune_df) > 0) {
                write.csv(results$local_qc_passed_exposure_variants_before_ld_prune_df, glue('{summary_folder}{exposure_id}_qc_passed_exposure_variants_before_ld_prune.csv'), row.names = FALSE)
            }
            if (nrow(results$local_qc_passed_exposure_variants_after_ld_prune_df) > 0) {
                write.csv(results$local_qc_passed_exposure_variants_after_ld_prune_df, glue('{summary_folder}{exposure_id}_qc_passed_exposure_variants_after_ld_prune.csv'), row.names = FALSE)
            }
            if (nrow(results$local_qc_passed_exposure_variants_failed_during_ld_prune_df) > 0) {
                write.csv(results$local_qc_passed_exposure_variants_failed_during_ld_prune_df, glue('{summary_folder}{exposure_id}_qc_passed_exposure_variants_failed_during_ld_prune.csv'), row.names = FALSE)
            }
            local_exposure_with_no_significant_variants_before_ld_prune <- c(local_exposure_with_no_significant_variants_before_ld_prune, results$local_exposure_with_no_significant_variants_before_ld_prune)
            local_exposure_with_no_significant_variants_after_ld_prune <- c(local_exposure_with_no_significant_variants_after_ld_prune, results$local_exposure_with_no_significant_variants_after_ld_prune)


            if (nrow(exposure_df) > 0) {
                for (outcome_id in outcome_file_df$gwas_id) {
                    local_outcome_with_no_variant_df <- data.frame(outcome_id = character(), exposure_id = character(), stringsAsFactors = FALSE)
                    dat <- data.frame()
                    local_harmonised_df <- data.frame()
                    results <- list()
                    outcomevcf <- glue("{output_folder}/filtered_outcome_vcf/{outcome_id}_filt.vcf.gz")
                    prefix <- glue("{out_folder}{exposure_id}_{outcome_id}")
                    results <- Harmonise_exposure_Outcome_data(outcome_id, exposure_id, exposure_df, outcomevcf, 
                                                               ld_file, pvalue_cutoff, tag_r2,
                                                               local_outcome_with_no_variant_df,
                                                               local_harmonised_df)
                    dat <- results$dat
                    if (nrow(results$local_outcome_with_no_variant_df) > 0) {
                        write.csv(results$local_outcome_with_no_variant_df, glue('{summary_folder}{exposure_id}_{outcome_id}_outcome_with_no_variant.csv'), row.names = FALSE)
                    }


                    if (nrow(dat) > 0) {
                        dat$`samplesize.outcome` <- 33000

                        mr_analysis_result <- list()
                        mr_single_variant_result <- list()
                        mr_res <- data.frame()
                        mendelian_randomization_df <- data.frame()
                        mendelian_randomization_ivw_delta_df <- data.frame()

                        # Two Sample MR test
                        mr_analysis_result <- perform_two_sample_mr_analysis(dat, outcome_id, exposure_id)
                        mr_single_variant_result <- perform_single_variant_mr_analysis(dat, outcome_id, exposure_id)

                        # Biogen mrpipeline
                        mr_res <- do_mr(dat, f_cutoff = 0, all_wr = TRUE, verbose = TRUE)

                        mendelian_randomization_df <- perform_mendelian_randomization(dat, outcome_id, exposure_id)
                        mendelian_randomization_ivw_delta_df <- perform_mendelian_randomization_ivw_delta(dat, outcome_id, exposure_id)

                        # Write CSV files
                        write.csv(results$local_harmonised_df, glue('{prefix}_Harmonised_Exposure_Outcome.csv'), row.names = FALSE)
                        write.csv(mr_analysis_result$mr_res, glue('{prefix}_TwoSampleMR_Analysis_Multiple_MR_Test.csv'), row.names = FALSE)
                        write.csv(mr_analysis_result$mr_het, glue('{prefix}_TwoSampleMR_Analysis_Heterogeneity_Test.csv'), row.names = FALSE)
                        write.csv(mr_analysis_result$mr_pleo, glue('{prefix}_TwoSampleMR_Analysis_Hpleiotropy_Test.csv'), row.names = FALSE)
                        write.csv(mr_analysis_result$mr_direction, glue('{prefix}_TwoSampleMR_Analysis_Directionality_Test.csv'), row.names = FALSE)
                        write.csv(mr_single_variant_result$mr_single_meta_fixed, glue('{prefix}_TwoSampleMR_Analysis_SingleVariantMetafixed_Test.csv'), row.names = FALSE)
                        write.csv(mr_single_variant_result$mr_single_wald_ratio, glue('{prefix}_TwoSampleMR_Analysis_SingleVariantWald_Test.csv'), row.names = FALSE)
                        write.csv(mr_res, glue('{prefix}_MendelianPipeline_Test.csv'), row.names = FALSE)
                        write.csv(mendelian_randomization_df, glue('{prefix}_MendelianRandomization_AllTest.csv'), row.names = FALSE)
                        write.csv(mendelian_randomization_ivw_delta_df, glue('{prefix}_MendelianRandomization_IVW_Delta_Test.csv'), row.names = FALSE)
                    }
                }


                for (suffix in  file_suffixes) {
                    m_df <- data.frame()
                    for (outcome_id in outcome_file_df$gwas_id) {
                        file_path <- glue("{out_folder}{exposure_id}_{outcome_id}_{suffix}.csv")
                        # Attempt to read the CSV file and concatenate to m_df
                        tryCatch({
                            if (file.exists(file_path)) {
                                df <- read.csv(file_path, sep = ",")
                                m_df <- bind_rows(m_df, df)
                            } else {
                                cat("File does not exist:", file_path, "\n")
                            }
                        }, error = function(e) { cat("Error reading or concatenating data:", conditionMessage(e), "\n") })
                    }
                    # Write the merged data frame to a CSV file
                    tryCatch({
                    if (suffix != "Harmonised_Exposure_Outcome") {
                        write.csv(m_df, glue('{out_folder}{exposure_id}_with_all_outcomes_{suffix}.csv'), row.names = FALSE)
                    } else {
                        write.csv(m_df, glue('{summary_folder}{exposure_id}_with_all_outcomes_{suffix}.csv'), row.names = FALSE)
                    }
                    }, error = function(e) { 
                    cat("Error writing merged data frame to CSV:", conditionMessage(e), "\n") 
                    })
                }
            } 
        } else {
                local_exposure_with_no_significant_variants_before_ld_prune <- c(local_exposure_with_no_significant_variants_before_ld_prune, exposure_id)
        }
            if (length(local_exposure_with_no_significant_variants_before_ld_prune) > 0) {
                write.csv(as.data.frame(local_exposure_with_no_significant_variants_before_ld_prune), glue('{summary_folder}{exposure_id}_exposure_with_no_significant_variants_before_ld_prune.csv'), row.names = FALSE)
            }
            if (length(local_exposure_with_no_significant_variants_after_ld_prune) > 0) {
                write.csv(as.data.frame(local_exposure_with_no_significant_variants_after_ld_prune), glue('{summary_folder}{exposure_id}_exposure_with_no_significant_variants_after_ld_prune.csv'), row.names = FALSE)
            }
    }

}


run_parallel_mendelian_analysis_pipeline <- function(exposure_file_df, outcome_file_df, output_folder, dbsnp_df, ld_file, 
                                                     infoscore_cutoff, maf_cutoff, pvalue_cutoff, variant_id_type, 
                                                     mhc_remove, mhc_region, remove_chrX, 
                                                     outcomevcf_path, tag_r2, n_vcf_parallel, N) {
    # Convert to data.table if it's not already one
    setDT(exposure_file_df)
    # Generate a group index for splitting
    exposure_file_df[, group := cut(seq_len(.N), N, labels = FALSE)]
    # Split the data.table into a list of data.tables based on the group
    job_batches <- split(exposure_file_df, by = "group")
    # Remove the group column from each split data.table
    job_batches <- lapply(job_batches, function(dt) dt[, group := NULL])
    
    # Parallel execution using mclapply exposure vcf files
    mclapply(1:length(job_batches), function(i) {
        mendelian_analysis_pipeline(job_batches[[i]], i, outcome_file_df, output_folder, dbsnp_df, ld_file, 
                                        infoscore_cutoff, maf_cutoff, pvalue_cutoff, variant_id_type, mhc_remove, 
                                        mhc_region, remove_chrX, outcomevcf_path, tag_r2)
    }, mc.cores = n_vcf_parallel)
    
        patterns <- c("exposure_with_no_significant_variants_before_ld_prune","exposure_with_no_significant_variants_after_ld_prune",
                    "qc_passed_exposure_variants_before_ld_prune","qc_passed_exposure_variants_after_ld_prune",
                    "outcome_with_no_variant","Harmonised_Exposure_Outcome")

        # Function to read and merge CSV files for each pattern
        summary_folder<-glue("{output_folder}/summary_exposure_outcome/")
        merge_and_save_csv <- function(pattern) {
            file_paths <- list.files(path = glue("{output_folder}/summary_exposure_outcome/"),
                                    pattern = paste0(pattern, "\\.csv$"), full.names = TRUE)
            if (length(file_paths) > 0) {
                merged_df <- bind_rows(lapply(file_paths, read.csv))
                output_file <- glue("{output_folder}/{pattern}_merged.csv")
                write.csv(merged_df, file = output_file, row.names = FALSE)
            } else {
                message(glue("No files found for pattern: {pattern}"))
            }
        }
        # Apply the function to each pattern
        lapply(patterns, merge_and_save_csv)
    }



# Parallel execution using mclapply outcome vcf files
results <- mclapply(1:nrow(outcome_file_df), function(i) {
    vcf_filter(i, outcome_file_df, lower_AF, upper_AF, infoscore_cutoff, "filtered_outcome_vcf", variant_id_type = variant_id_type, p_value = 0.000001,variants_include=variants_include)
}, mc.cores = n_vcf_parallel)


#Reading dbsnp file 
tryCatch({
    if (dbsnpfile != "NA") {
        dbsnp_df <- fread(dbsnpfile)
        colnames(dbsnp_df) <- c('seqnames', 'start', 'ID', 'REF', 'ALT')
        dbsnp_df$seqnames <- as.numeric(dbsnp_df$seqnames)
        dbsnp_df$seqnames <- as.factor(dbsnp_df$seqnames)
    } else {
        dbsnp_df <- "NA"
    }
}, error = function(e) {
    stop(glue("Error processing dbSNP file: {conditionMessage(e)}"))
})


# Parallel execution using mclapply exposure vcf files
results <- mclapply(1:nrow(exposure_file_df), function(i) {
    vcf_filter(i, exposure_file_df, lower_AF, upper_AF, infoscore_cutoff, "filtered_exposure_vcf",
     variant_id_type = variant_id_type, p_value = -log10(pvalue_cutoff),variants_include=variants_include)
}, mc.cores = n_vcf_parallel)


# Parallel execution using mclapply run_parallel_mendelian_analysis_pipeline 
run_parallel_mendelian_analysis_pipeline(exposure_file_df, outcome_file_df, output_folder, dbsnp_df, ld_file, 
                                          infoscore_cutoff, maf_cutoff, pvalue_cutoff, variant_id_type, mhc_remove, 
                                          mhc_region, remove_chrX, outcomevcf_path, tag_r2, n_vcf_parallel, N)
