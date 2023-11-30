setwd("/home/mdiaz/lungsquamouscells/prepo_data")
library(TCGAbiolinks)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)
library(dplyr)

# Load the TxDb annotation database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Get information about promoters
promoters <- promoters(txdb, upstream = 1000, downstream = 0)

# Visualize the results
head(promoters)

# Path to the BED file of enhancers downloaded from ENCODE
enhancer_file <- "/home/mdiaz/Descargas/ENCFF596CUU.bed.gz"

# Load enhancers from the BED file
enhancers <- fread(enhancer_file, header = FALSE)

# Select only the first three columns
enhancers <- enhancers[, 1:3, with = FALSE]

# Assign names to the columns
colnames(enhancers) <- c("chr", "start", "end")

# Visualize the first results
head(enhancers)

##################################
# Load enhancers from the BED file excluding the header
enhancers <- fread(enhancer_file, skip = 1)

# Select only the first three columns
enhancers <- enhancers[, 1:3, with = FALSE]

# Assign names to the columns
colnames(enhancers) <- c("chr", "start", "end")

# Visualize the first results
head(enhancers)

##########################
# Convert to GRanges formats for intersection use
promoters_gr <- makeGRangesFromDataFrame(promoters, start.field = "start", end.field = "end")
enhancers_gr <- makeGRangesFromDataFrame(enhancers, start.field = "start", end.field = "end")

# Print the lengths of the datasets
cat("Number of promoters:", length(promoters_gr), "\n")
cat("Number of enhancers:", length(enhancers_gr), "\n")

# Create data frames
promoters_df <- as.data.frame(promoters_gr)
enhancers_df <- as.data.frame(enhancers_gr)

# Add information about the element type (promoter or enhancer)
promoters_df$type <- "promoter"
enhancers_df$type <- "enhancer"

# Combine data frames
combined_df <- rbind(promoters_df, enhancers_df)

# Create data frames 
promoters_df <- as.data.frame(promoters_gr)
enhancers_df <- as.data.frame(enhancers_gr)

# Add information about the element type (promoter or enhancer)
promoters_df$type <- "promoter"
enhancers_df$type <- "enhancer"

# Combine data frames
combined_df <- rbind(promoters_df, enhancers_df)

combined_df$ENST_ID <- rownames(combined_df)

promo_enhan <- c("ENST_ID", "seqnames", "start", "end", "width", "strand", "type")

selec_promoters_enhancers <- combined_df[promo_enhan]
selec_promoters_enhancers2 <- selec_promoters_enhancers
selec_promoters_enhancers2$ENST_ID <- sub("\\.\\d+$", "", selec_promoters_enhancers2$ENST_ID)

# Print the updated data frame
print(selec_promoters_enhancers2)

# Save the dataframe as a TSV file
write.table(selec_promoters_enhancers, "promoters_enhancers.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


###########################
# ENSEMBL IDs for CpG sites
###########################

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# List of genes of interest
genes <- c("ZNF334", "FLT3", "KCNB1", "ARHGDIG", "H2BC3", "C8orf76", "SCGB3A1", "GPR27",
           "ZNF542P", "GAL3ST2", "BCOR", "TBXAS1", "SRCIN1", "TCEA3", "GMDS", "ABP3",
           "MYADML", "PNMA8A", "COL6A6", "ZFHX4", "PLPPR3", "KRT7", "ACTB", "SFMBT1")

# Get the ENST IDs for each gene
get_enst <- function(gene_symbol) {
  gene_info <- getBM(attributes = c("ensembl_transcript_id", "external_gene_name"),
                     filters = "external_gene_name",
                     values = gene_symbol,
                     mart = ensembl)
  return(gene_info)
}

# Create a data frame with the results
result_df <- data.frame()
for (gene in genes) {
  gene_info <- get_enst(gene)
  result_df <- rbind(result_df, gene_info)
}

# Print the result
print(result_df)

colnames(result_df) <- c("ENST_ID", "external_gene_name")
#########################################################
#### MERGE CpG SITES WITH PROMOTERS & ENHANCERS DATABASE#
#########################################################

ensembl_CpGsites <- result_df %>% 
  left_join(selec_promoters_enhancers2, by=c("ENST_ID"))

ensembl_CpGsites$ENST <- NULL

colnames(ensembl_CpGsites) <- c("ENST_ID", "external_CpG_name", "seqnames", "start", "end", "width", "strand", "type")

write.table(ensembl_CpGsites, "ensembl_CpGsites.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
