#!/usr/bin/env Rscript

library("DESeq2")
library(tidyverse)
library(plyr)
#install.packages("limma")
library(limma)

reading_gene_counts_folder <- function(input_folder, outputFolder) {
  file.name=c(list.files(input_folder, pattern = "gene$", full.names = T))
  for (i in 1:length(file.name)){
    table.file=read_tsv(file.name[i], skip = 1)
    Count= colnames(table.file)[7]
    #table.file=table.file[1:10,]
    if (i==1){
      counts=as.data.frame(table.file[,c(1,2,3,4,5,6,7)])
      vec=as.character(c("Geneid", "Chr", "Start", "End", "Strand", "Length", Count))}
    if (i>1){
      counts=cbind(counts,table.file[,7])
      vec=c(vec,Count)
      vec=gsub("_STAR_aln_to_genome_q10_srt.bam", "", vec)
    }}
  colnames(counts)=vec
  
  #Import gene name conversion table used for 3prime RNA seq
  read_tsv("/Users/chill/Library/CloudStorage/OneDrive-TheWistarInstitute/professional/git_hub/scripts_to_upload/gsva/Homo_sapiens.GRCh37.87.conversion.txt", col_names = F) %>% 
    dplyr::rename(chr=X1, start=X2, end = X3, name = X4, transcript_id = X5, strand = X6, gene_id = X7) -> ens_hgnc
  
  #View table of genes with name conversions
  inner_join(counts, ens_hgnc, by = c("Geneid"="gene_id")) %>% 
    dplyr::select(-c("chr","start","end", "transcript_id", "strand")) %>% 
    dplyr::select(name, 1:ncol(.)) %>% 
    dplyr::rename(Gene_name="name") %>% 
    base::unique()
  
  inner_join(counts, ens_hgnc, by = c("Geneid"="gene_id")) %>% 
    select(-c("chr","start","end", "transcript_id", "strand")) %>% 
    select(name, 1:ncol(.)) %>% 
    dplyr::rename(Gene_name="name") %>% 
    base::unique() -> gene_name_counts
  
  #Make DeSeq nreads table 
  counts %>% 
    #Omit any genes with zero counts across samples:
    mutate(sum=rowSums(.[c(-1,-2,-3,-4,-5,-6)])) %>%
    filter(sum > 0) %>% 
    dplyr::select(-c(2:6),-sum) %>%
    column_to_rownames('Geneid') %>%
    +1 -> Nreads
  
  tpm <- function(counts,len) {
    x <- counts/len
    return(t(t(x)*1e6/colSums(x)))
  }
  
  g <- data.frame( ensembl_gene_id = counts$Geneid , 
                   transcript_length = counts$Length,
                   stringsAsFactors = F, row.names = counts$Geneid)
  
  g <- g[!duplicated(g$ensembl_gene_id),]
  
  igenes <- intersect(rownames(Nreads),g$ensembl_gene_id)
  g1 <- g[igenes,]
  cf1 <- Nreads[igenes,]
  all.equal(rownames(cf1),g1$ensembl_gene_id)
  
  ct <- tpm(cf1,g1$transcript_length)
  ct <- log2( ct + 1 )
  #boxplot(ct,ylab=expression('Log'[2]~'Read counts'),las=2,main="TPM")
  
  #this is from https://nbisweden.github.io/workshop-RNAseq/2011/lab_preprocessing.html#31_CPMTPM
  
  file_name <- "counts_tpm.txt"
  
  ct %>% data.frame() %>% rownames_to_column(var = "gene_id") -> ct
  
  write.table(ct, file = file.path(outputFolder, file_name), col.names = TRUE, row.names = FALSE, sep = "\t")
  
  

}

reading_gene_counts_files <- function(input_file_list, outputFolder) {
  file.name=c(input_file_list)
  for (i in 1:length(file.name)){
    table.file=read_tsv(file.name[i], skip = 1)
    Count= colnames(table.file)[7]
    #table.file=table.file[1:10,]
    if (i==1){
      counts=as.data.frame(table.file[,c(1,2,3,4,5,6,7)])
      vec=as.character(c("Geneid", "Chr", "Start", "End", "Strand", "Length", Count))}
    if (i>1){
      counts=cbind(counts,table.file[,7])
      vec=c(vec,Count)
      vec=gsub("_STAR_aln_to_genome_q10_srt.bam", "", vec)
    }}
  colnames(counts)=vec
  
  #Import gene name conversion table used for 3prime RNA seq
  read_tsv("/Users/chill/Library/CloudStorage/OneDrive-TheWistarInstitute/professional/git_hub/scripts_to_upload/gsva/Homo_sapiens.GRCh37.87.conversion.txt", col_names = F) %>% 
    dplyr::rename(chr=X1, start=X2, end = X3, name = X4, transcript_id = X5, strand = X6, gene_id = X7) -> ens_hgnc
  
  #View table of genes with name conversions
  inner_join(counts, ens_hgnc, by = c("Geneid"="gene_id")) %>% 
    dplyr::select(-c("chr","start","end", "transcript_id", "strand")) %>% 
    dplyr::select(name, 1:ncol(.)) %>% 
    dplyr::rename(Gene_name="name") %>% 
    base::unique()
  
  inner_join(counts, ens_hgnc, by = c("Geneid"="gene_id")) %>% 
    select(-c("chr","start","end", "transcript_id", "strand")) %>% 
    select(name, 1:ncol(.)) %>% 
    dplyr::rename(Gene_name="name") %>% 
    base::unique() -> gene_name_counts
  
  #Make DeSeq nreads table 
  counts %>% 
    #Omit any genes with zero counts across samples:
    mutate(sum=rowSums(.[c(-1,-2,-3,-4,-5,-6)])) %>%
    filter(sum > 0) %>% 
    dplyr::select(-c(2:6),-sum) %>%
    column_to_rownames('Geneid') %>%
    +1 -> Nreads
  
  tpm <- function(counts,len) {
    x <- counts/len
    return(t(t(x)*1e6/colSums(x)))
  }
  
  g <- data.frame( ensembl_gene_id = counts$Geneid , 
                   transcript_length = counts$Length,
                   stringsAsFactors = F, row.names = counts$Geneid)
  
  g <- g[!duplicated(g$ensembl_gene_id),]
  
  igenes <- intersect(rownames(Nreads),g$ensembl_gene_id)
  g1 <- g[igenes,]
  cf1 <- Nreads[igenes,]
  all.equal(rownames(cf1),g1$ensembl_gene_id)
  
  ct <- tpm(cf1,g1$transcript_length)
  ct <- log2( ct + 1 )
  #boxplot(ct,ylab=expression('Log'[2]~'Read counts'),las=2,main="TPM")
  
  #this is from https://nbisweden.github.io/workshop-RNAseq/2011/lab_preprocessing.html#31_CPMTPM
  
  file_name <- "counts_tpm.txt"
  
  ct %>% data.frame() %>% rownames_to_column(var = "gene_id") -> ct
  
  write.table(ct, file = file.path(outputFolder, file_name), col.names = TRUE, row.names = FALSE, sep = "\t")
  
}

#to stop parsing when it encounters another --Flag 
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

#print(args)
# Process --outputFolder argument
outputFolder_args <- parseArgs(args, "--outputFolder")

print(outputFolder_args)

# Process --QuantseqGeneFolder argument
Quantseq_tpm_folder_args <- parseArgs(args, "--QuantseqGeneFolder")
if (length(Quantseq_tpm_folder_args) >= 1) {
  reading_gene_counts_folder(Quantseq_tpm_folder_args, outputFolder_args)
}

# Process --QuantseqSampleTpmCounts argument
Quantseq_tpm_args <- parseArgs(args, "--QuantseqGeneFiles")
if (length(Quantseq_tpm_args) > 0) {
  reading_gene_counts_files(Quantseq_tpm_args, outputFolder_args)
}







