#!/usr/bin/env Rscript

library(tidyverse)
library(GSVA)
library(dplyr)
library(gplots)
library(reshape2)
library(msigdbr)


############################################################################################################this portion of the script is if you want the  msgidbr collection in your gene signatures to test

msigdbrCollectionAdding <- function(speciesName, collectionName) {
  #msigdbr_species() list species if needed 
  hs_gsea <- msigdbr(species = speciesName) #download specicies specific collection
  hs_gsea %>% 
    dplyr::distinct(gs_cat, gs_subcat) %>% 
    dplyr::arrange(gs_cat, gs_subcat)
  hs_gsea_c2 <- msigdbr(species = speciesName, # change depending on species your data came from
                        category = collectionName) %>% # choose your msigdb collection of interest
    dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in
  return(hs_gsea_c2)
  
}




############################################################################################################
#this portion of the script is if you want a folders your own collections in your gene signatures to test

readingInFolder <- function(directory) {
  # Get the list of files ending in _nearest_gene.bed
  file_list <- list.files(directory, pattern = "_nearest_gene\\.bed$", full.names = TRUE)
  
  # Create an empty list to store the modified data frames
  modified_data_list <- list()
  
  # Iterate over each file
  for (file in file_list) {
    # Check if the file is empty
    if (file.size(file) == 0) {
      next  # Skip to the next file if it's empty
    }
    
    # Read the file
    data <- read.table(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    
    # Randomly select 500 rows (or fewer if there are fewer than 500 rows)
    num_rows <- min(nrow(data), 500)
    selected_rows <- sample(1:nrow(data), size = num_rows, replace = FALSE)
    selected_data <- data[selected_rows, ]
    
    # Extract the gene symbol from the 3rd to last column
    gene_symbol <- selected_data[, ncol(selected_data) - 2]
    
    # Create a new data frame with only the gene_symbol and gs_name columns
    modified_data <- data.frame(gs_name = basename(file), gene_symbol = gene_symbol)
    
    # Store the modified data frame in the list
    modified_data_list[[basename(file)]] <- modified_data
  }
  
  if (length(modified_data_list) > 0) {
    # Combine all modified data frames into a single data frame
    combined_data_from_folders_C2 <- do.call(rbind, modified_data_list)
    
    # Remove row names from the combined data frame
    rownames(combined_data_from_folders_C2) <- NULL
    
    # Return the combined data frame
    return(combined_data_from_folders_C2)
  } else {
    # Return NULL if no valid files found
    return(NULL)
  }
}



############################################################################################################
#this portion of the script is if you want add specific collections in your gene signatures to test

readingInFiles <- function(file_list) {
  # Create an empty list to store the modified data frames
  modified_data_list <- list()
  
  # Iterate over each file in the provided list
  for (file in file_list) {
    # Check if the file is empty
    if (file.size(file) == 0) {
      next  # Skip to the next file if it's empty
    }
    
    # Read the file
    data <- read.table(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    
    # Randomly select 500 rows (or fewer if there are fewer than 500 rows)
    num_rows <- min(nrow(data), 500)
    selected_rows <- sample(1:nrow(data), size = num_rows, replace = FALSE)
    selected_data <- data[selected_rows, ]
    
    # Extract the gene symbol from the 3rd to last column
    gene_symbol <- selected_data[, ncol(selected_data) - 2]
    
    # Create a new data frame with only the gene_symbol and gs_name columns
    modified_data <- data.frame(gs_name = basename(file), gene_symbol = gene_symbol)
    
    # Store the modified data frame in the list
    modified_data_list[[basename(file)]] <- modified_data
  }
  
  if (length(modified_data_list) > 0) {
    # Combine all modified data frames into a single data frame
    combined_data_from_files_C2 <- do.call(rbind, modified_data_list)
    
    # Remove row names from the combined data frame
    rownames(combined_data_from_files_C2) <- NULL
    
    # Return the combined data frame
    return(combined_data_from_files_C2)
  } else {
    # Return NULL if no valid files found
    return(NULL)
  }
}

############################################################################################################
#this portion of the script is if you want use TPM gene expression counts downloaded form the cBIOportal and combined into a matrix using expression_matrix_creation_from_cbioportal.py in your sample counts to test
readingInCbioCounts <- function(combined_gene_countsFile) {
  
  
  data <- read.table(combined_gene_countsFile, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  
  colnames(data) <- data[1, ] 
  
  data <- data[- 1, ]
  
  
  # Remove duplicate rows based on gene_name
  data <- data[!duplicated(data$gene_name), ]
  
  # Exclude specific columns from the data frame
  columns_to_exclude <- c("gene_id", "gene_type")
  CbioCounts_TPMdf <- data[, !(names(data) %in% columns_to_exclude), drop = FALSE]
  
  return(CbioCounts_TPMdf)
}

############################################################################################################
#this portion of the script is if you want use TPM gene expression form the you own RNA-seq samples with Gene_ID's like ENSG00000240361.1 in first column and combined into a matrix using  in your sample counts to test
readingInRNAseqTPMCounts <- function(RNAseqTPMCountsFileList) {
  # Create an empty list to store the modified data frames
  modified_data_list <- list()
  
  #read in list of genes and gene ids 
  read_tsv("/Users/chill/Library/CloudStorage/OneDrive-TheWistarInstitute/professional/git_hub/scripts_to_upload/gsva/gencodeV24lift37_basicbed.txt", col_names = F) %>% 
    dplyr::rename(chr=X1, start=X2, end = X3, name = X4, transcript_id = X5, strand = X6, gene_id = X7) -> ens_hgnc_regular
  
  # Iterate over each file in the provided list
  for (file in RNAseqTPMCountsFileList) {
    # Check if the file is empty
    if (file.size(file) == 0) {
      next  # Skip to the next file if it's empty
    }
    
    # Read the file
    data <- read.table(file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    # Store the modified data frame in the list
    modified_data_list[[basename(file)]] <- data
  }
  
  if (length(modified_data_list) > 0) {
    # Combine all modified data frames into a single data frame
    
    merged_RNAseq_TPMdf <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), modified_data_list)
    
    merge(merged_RNAseq_TPMdf, ens_hgnc_regular) -> merged_RNAseq_TPMdf
    
    #print(merged_RNAseq_TPMdf)
    
    merged_RNAseq_TPMdf <- merged_RNAseq_TPMdf %>% dplyr::select(-start,-end, -transcript_id, -strand, -chr, -gene_id)
    
    merged_RNAseq_TPMdf <- merged_RNAseq_TPMdf[!duplicated(merged_RNAseq_TPMdf$name), ]
    
    merged_RNAseq_TPMdf %>% dplyr::rename(gene_name = name) -> merged_RNAseq_TPMdf
    
    
    # Return the combined data frame
    return(merged_RNAseq_TPMdf)
  } else {
    # Return NULL if no valid files found
    return(NULL)
  }
}

############################################################################################################
#this portion of the script is if you want use TPM gene expression form the you own RNA-seq samples with Gene_ID's like ENSG00000240361.1 in first column and combined into a matrix using  in your sample counts to test
readingInQuantseqTPMCounts <- function(QuantseqTPMCountsFileList) {
  # Create an empty list to store the modified data frames
  modified_data_list <- list()
  
  #read in list of genes and gene ids 
  read_tsv("/Users/chill/Library/CloudStorage/OneDrive-TheWistarInstitute/professional/git_hub/scripts_to_upload/gsva/Homo_sapiens.GRCh37.87.conversion.txt", col_names = F) %>% 
    dplyr::rename(chr=X1, start=X2, end = X3, name = X4, transcript_id = X5, strand = X6, gene_id = X7) -> ens_hgnc_regular
  
  # Iterate over each file in the provided list
  for (file in QuantseqTPMCountsFileList) {
    # Check if the file is empty
    if (file.size(file) == 0) {
      next  # Skip to the next file if it's empty
    }
    
    # Read the file
    data <- read.table(file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    # Store the modified data frame in the list
    modified_data_list[[basename(file)]] <- data
  }
  
  if (length(modified_data_list) > 0) {
    # Combine all modified data frames into a single data frame
    
    merged_Quantseq_TPMdf <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), modified_data_list)
    
    
    merge(merged_Quantseq_TPMdf, ens_hgnc_regular) -> merged_Quantseq_TPMdf
    
    merged_Quantseq_TPMdf <- merged_Quantseq_TPMdf %>% dplyr::select(-start,-end, -transcript_id, -strand, -chr, -gene_id)
    
    merged_Quantseq_TPMdf <- merged_Quantseq_TPMdf[!duplicated(merged_Quantseq_TPMdf$name), ]
    
    merged_Quantseq_TPMdf %>% dplyr::rename(gene_name = name) -> merged_Quantseq_TPMdf
    
    
    # Return the combined data frame
    return(merged_Quantseq_TPMdf)
  } else {
    # Return NULL if no valid files found
    return(NULL)
  }
}


#too summarize data with standard devations
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#too plot and filter by primary diagnosis column in the clinicalDataSet
filter_by_set_names <- function(set_names) {
  for (selected_set_name in set_names) {
    merged_df %>% filter(set_name == selected_set_name) -> merged_df_set_name_filtered
    merged_df_set_name_filtered <- distinct(merged_df_set_name_filtered, case_id, .keep_all = TRUE)
    df2 <- data_summary(merged_df_set_name_filtered, varname="gsva_score", groupnames=c("primary_diagnosis"))
    df2$up_down <- ifelse(df2$gsva_score < 0, "blue", "red")
    file_name <- paste(selected_set_name, "_filtered_merged_summarized_by_set_name_primary_diagnosis_results.txt", sep = "_")
    write.table(df2, file = file.path(outputFolder, file_name), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    file_name <- paste(selected_set_name, "_filtered_merged_summarized_by_set_name_primary_diagnosis_results.pdf", sep = "_")
    df2 %>% slice_max(order_by = gsva_score, n = 10) -> df_top_10
    df2 %>% slice_min(order_by = gsva_score, n = 10) -> df_bottom_10
    df_2_top_bottom_10 <- bind_rows(df_top_10, df_bottom_10)
    ggplot(df_2_top_bottom_10, aes(x = reorder(primary_diagnosis, gsva_score), y = gsva_score)) +
      geom_bar(position = "dodge", stat = "identity", fill = df_2_top_bottom_10$up_down) +
      geom_errorbar(aes(x = reorder(primary_diagnosis, gsva_score), ymin = gsva_score - sd, ymax = gsva_score + sd),
                    width = 0.2, position = position_dodge(0.9)) +
      coord_flip() +
      theme_classic() +
      xlab("primary_diagnosis") +
      ylab(gsvaMethod) +
      theme(legend.position = "bottom")
    
    ggsave(filename = file.path(outputFolder, file_name), width = 10, height = 10)
  }
}


#too plot and filter by gene set
filter_by_primary_diagnosis <- function(set_primary_diagnosis) {
  for (selected_primary_diagnosis in set_primary_diagnosis) {
    merged_df %>% filter(primary_diagnosis == "cell_line1") -> merged_df_set_name_filtered
    df2 <- data_summary(merged_df_set_name_filtered, varname="gsva_score", groupnames=c("set_name"))
    df2$up_down <- ifelse(df2$gsva_score < 0, "blue", "red")
    file_name <- paste(selected_primary_diagnosis, "_filtered_merged_summarized_by_set_name_primary_diagnosis_results.txt", sep = "_")
    write.table(df2, file = file.path(outputFolder, file_name), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    file_name <- paste(selected_primary_diagnosis, "_filtered_merged_summarized_by_set_name_primary_diagnosis_results.pdf", sep = "_")
    df2 %>% slice_max(order_by = gsva_score, n = 10) -> df_top_10
    df2 %>% slice_min(order_by = gsva_score, n = 10) -> df_bottom_10
    df_2_top_bottom_10 <- bind_rows(df_top_10, df_bottom_10)
    ggplot(df_2_top_bottom_10, aes(x = reorder(set_name, gsva_score), y = gsva_score)) +
      geom_bar(position = "dodge", stat = "identity", fill = df_2_top_bottom_10$up_down) +
      geom_errorbar(aes(x = reorder(set_name, gsva_score), ymin = gsva_score - sd, ymax = gsva_score + sd),
                    width = 0.2, position = position_dodge(0.9)) +
      coord_flip() +
      theme_classic() +
      xlab("set_name") +
      ylab(gsvaMethod) +
      theme(legend.position = "bottom")
    ggsave(filename = file.path(outputFolder, file_name), width = 10, height = 10)
  }
}

isFlag <- function(x) {
  grepl("^--", x)
}


# Helper function to parse arguments
parseArgs <- function(args, flag) {
  flag_index <- which(args == flag)
  if (length(flag_index) > 0 && flag_index < length(args)) {
    next_args <- args[(flag_index + 1):length(args)]
    next_flag_index <- which(isFlag(next_args))
    if (length(next_flag_index) > 0) {
      return(next_args[1:(next_flag_index - 1)])
    } else {
      return(next_args)
    }
  }
  return(character(0)) # Return empty character vector if flag not found or no values after flag
}

# Initialize variables
gsvaMethod <- NULL

# Get all command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Process --analysisMethod argument
analysis_method_args <- parseArgs(args, "--analysisMethod")
if (length(analysis_method_args) >= 1) {
  gsvaMethod <- analysis_method_args[1]
  valid_methods <- c("ssgsea", "gsva", "zscore", "plage")
  if (!gsvaMethod %in% valid_methods) {
    cat("Wrong value for analysisMethod. Accepted values are: ssgsea, gsva, zscore, or plage\n")
    quit(status = 1) # Exit the script if invalid method is provided
  }
}



# Process --msigdbr argument
msigdbr_args <- parseArgs(args, "--msigdbr")
if (length(msigdbr_args) >= 2) {
  hs_gsea_c2 <- msigdbrCollectionAdding(msigdbr_args[1], msigdbr_args[2])
}


# Process --folderNearestGeneBed argument
folderNearestGeneBed_args <- parseArgs(args, "--folderNearestGeneBed")
if (length(folderNearestGeneBed_args) >= 1) {
  combined_data_from_folders_c2 <- readingInFolder(folderNearestGeneBed_args[1])
}

# Process --specificNearestGeneBed argument
specific_gene_args <- parseArgs(args, "--specificNearestGeneBed")
#cat("specific_gene_args results: ", specific_gene_args)
if (length(specific_gene_args) > 0) {
  combined_data_from_files_c2 <- readingInFiles(specific_gene_args)
}


# Process --cBioCounts argument
cbio_counts_args <- parseArgs(args, "--cBioCounts")
#cat("cbio_counts_args results: ",cbio_counts_args)
if (length(cbio_counts_args) >= 1) {
  CbioCounts_TPMdf <- readingInCbioCounts(cbio_counts_args[1])
}


# Process --RNAseqSampleTpmCounts argument
rnaseq_tpm_args <- parseArgs(args, "--RNAseqSampleTpmCounts")
#cat("rnaseq_tpm_args results: ", rnaseq_tpm_args)
if (length(rnaseq_tpm_args) > 0) {
  merged_RNAseq_TPMdf <- readingInRNAseqTPMCounts(rnaseq_tpm_args)
}


# Process --QuantseqSampleTpmCounts argument
Quantseq_tpm_args <- parseArgs(args, "--QuantseqSampleTpmCounts")
#cat("rnaseq_tpm_args results: ", rnaseq_tpm_args)
if (length(Quantseq_tpm_args) > 0) {
  merged_Quantseq_TPMdf <- readingInQuantseqTPMCounts(Quantseq_tpm_args)
}

# Process --outputFolder argument
outputFolder_args <- parseArgs(args, "--outputFolder")
#cat("rnaseq_tpm_args results: ", rnaseq_tpm_args)
if (length(outputFolder_args) >= 1) {
  outputFolder <- outputFolder_args[1]
}


# Process --clinicalDataSet argument
clinicalDataSet_args <- parseArgs(args, "--clinicalDataSet")
#cat("rnaseq_tpm_args results: ", rnaseq_tpm_args)
if (length(clinicalDataSet_args) >= 1) {
  clinicalDataSet <- clinicalDataSet_args[1]
}






#going to combine gene sets 



c2_list <- mget(ls(pattern = "_c2"))

#print(c2_list)

dplyr::bind_rows(c2_list) -> hs_gsea_c2_combined



# Initialize an empty list
result_list <- list()

# Iterate over unique gene symbols
unique_symbols <- unique(hs_gsea_c2_combined$gs_name)
for (symbol in unique_symbols) {
  # Subset the data frame based on the gene symbol
  subset_df <- hs_gsea_c2_combined[hs_gsea_c2_combined$gs_name == symbol, ]
  
  # Get the gene_symbol values for the current gene symbol
  gs_names <- subset_df$gene_symbol
  
  # Add the gs_names to the result list with the gene symbol as the name
  result_list[[symbol]] <- gs_names
}

#lapply(result_list, head)

#now combining the TPM counts 

dfTPM_list <- mget(ls(pattern = "_TPMdf"))


all_TPM_combined <- Reduce(function(x, y) merge(x, y, by = "gene_name", all = TRUE), dfTPM_list)

all_TPM_combined  %>% drop_na() -> all_TPM_combined


file_name <- "all_TPM_combined.txt"
write.table(all_TPM_combined, file = file.path(outputFolder, file_name), sep = "\t", row.names = FALSE)


# Convert the data frame to a matrix and set the row names
gene_matrix <- as.matrix(all_TPM_combined[, -1])
rownames(gene_matrix) <- all_TPM_combined$gene_name

class(gene_matrix) <- "numeric"



#now doing the analysis 


#now running the gsva 

gbm_es <-  GSVA::gsva(gene_matrix, #your data
                      result_list, #signatures
                      min.sz=5, max.sz=500, #criteria for filtering gene sets
                      mx.diff=FALSE,
                      method="ssgsea") #options for method are "gsva", ssgsea', "zscore" or "plage"

#export matrix 


file_name <- paste(gsvaMethod, "results.txt", sep = "_")
write.table(gbm_es, file = file.path(outputFolder, file_name), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

#transform the matrix into data frame where one row is the gene_name (from the row names), another is the case_id (from the column names), the other column is named gsva_score (from the data)


df <- melt(gbm_es, id.vars = "gene_name", variable.name = "case_id", value.name = "gsva_score")

colnames(df) <- c("set_name", "case_id", "gsva_score")

# Read the clinical data into a data frame
clinical_data <- read.table(clinicalDataSet, header = TRUE, stringsAsFactors = FALSE, sep = "\t", , fill = TRUE, quote = "")


# Select specific columns from the clinical data
selected_columns <- c("case_id", "project_id", "primary_diagnosis", "site_of_resection_or_biopsy", "tissue_or_organ_of_origin")
clinical_selected <- clinical_data[selected_columns]

# Merge the clinical data with the existing data frame based on case_id
merged_df <- merge(df, clinical_selected, by = "case_id", all.x = TRUE)

# Print the resulting merged data frame
#print(merged_df)

merged_df$up_down <- ifelse(merged_df$gsva_score < 0, "Downregulated", "Upregulated")

file_name <- paste(gsvaMethod, "merged_results.txt", sep = "_")

write.table(merged_df, file = file.path(outputFolder, file_name), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(merged_df, varname="gsva_score", 
                    groupnames=c("set_name", "primary_diagnosis"))

df2$up_down <- ifelse(df2$gsva_score < 0, "blue", "red")

df2_filtered <- df2 %>%
  filter(!is.na(sd))

file_name <- paste(gsvaMethod, "merged_summarized_by_set_name_primary_diagnosis_results.txt", sep = "_")

write.table(df2_filtered, file = file.path(outputFolder, file_name), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



# Process --setNameToFilter argument
setNameToFilter_args <- parseArgs(args, "--setNameToFilter")
#cat("rnaseq_tpm_args results: ", rnaseq_tpm_args)
if (length(setNameToFilter_args) > 0) {
  filter_by_set_names(setNameToFilter_args)
}

# Process --setPrimaryDiagnosisFilter argument
primaryDiagnosisFilter_args <- parseArgs(args, "--primaryDiagnosisFilter")

print(primaryDiagnosisFilter_args)
#cat("rnaseq_tpm_args results: ", rnaseq_tpm_args)
if (length(primaryDiagnosisFilter_args) > 0) {
  filter_by_primary_diagnosis(primaryDiagnosisFilter_args)
}



