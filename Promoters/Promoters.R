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

#Count the frequency of each value
number_promoters_enhancers <- selec_promoters_enhancers2 %>%
  group_by(type) %>%
  summarise(count = n())

print(head(number_promoters_enhancers))

#Unique names for IDs
unique_ids <- unique(selec_promoters_enhancers2$ENST_ID)
length_ids <- length(unique_ids)
print(unique_ids)
print(length_ids)
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
            "MYADML", "PNMA8A", "COL6A6", "ZFHX4", "PLPPR3", "KRT7", "TRIO", "SFMBT1")
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

# Add the "main_gene_of_the_network" column
add_main_gene_column <- function(df) {
  df %>%
    mutate(
      network_main_gen = case_when(
        external_gene_name %in% c("ZNF334", "FLT3") ~ "mir-125b-2",
        external_gene_name %in% c("ARHGDIG", "KCNB1", "H2BC3", "C8orf76") ~ "mir-199a-2",
        external_gene_name %in% c("SCGB3A1", "GPR27", "GAL3ST2", "ZNF542P") ~ "CAPN2",
        external_gene_name %in% c("BCOR", "TBXAS1") ~ "LRP1",
        external_gene_name == "SRCIN1" ~ "RPS6",
        external_gene_name %in% c("TCEA3", "MYADML", "GMDS") ~ "RPS18",
        external_gene_name %in% c("PNMA8A", "COL6A6", "ZFHX4") ~ "EIF4G1",
        external_gene_name == "PLPPR3" ~ "EIF4G1/ACTB",
        external_gene_name == "KRT7" ~ "SDC1",
        external_gene_name == "TRIO" ~ "LRP1/DST",
        external_gene_name == "SFMBT1" ~ "DST",
        TRUE ~ NA_character_
      )
    )
}

result_df$main_gene_of_the_network <- NULL
# Add the "main_gene_of_the_network" column
result_df <- add_main_gene_column(result_df)

# Print the updated result
print(result_df)

result_df <- as.data.frame(result_df)
#########################################################
#### MERGE CpG SITES WITH PROMOTERS & ENHANCERS DATABASE#
#########################################################

ensembl_CpGsites <- result_df %>%
   left_join(selec_promoters_enhancers2, by=c("ENST_ID"))

ensembl_CpGsites$ENST <- NULL

colnames(ensembl_CpGsites) <- c("CpG_ENST_ID", "external_CpG_name", "network_main_gene", "CpG_seqnames", "CpG_start", "CpG_end", "CpG_width", "strand", "type")

write.table(ensembl_CpGsites, "ensembl_CpGsites.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

########################################
###ENSEMBL IDs FOR GENE NODES###########
########################################

#Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#List of genes of interest
genes <- c("CAPN2", "PFN2", "LRP1", "RPS6", "RPS18", "EIF4G1", "ACTB", "EEF2", "DST")

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
#### MERGE GENES WITH PROMOTERS & ENHANCERS DATABASE#####
#########################################################

ensembl_genes <- result_df %>%
  left_join(selec_promoters_enhancers2, by=c("ENST_ID"))

ensembl_genes$type <- NULL

colnames(ensembl_genes) <- c("ENST_ID", "external_gene_name", "seqnames", "start", "end", "width", "strand")

write.table(ensembl_genes, "ensembl_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


#####################################################
###Calculate distances between CAPN2 and CpG sites###
#####################################################
# Convert start and end positions to numeric values
ensembl_genes$start <- as.numeric(ensembl_genes$start)
ensembl_genes$end <- as.numeric(ensembl_genes$end)

ensembl_CpGsites$start <- as.numeric(ensembl_CpGsites$start)
ensembl_CpGsites$end <- as.numeric(ensembl_CpGsites$end)

# Function to calculate the midpoint for a given gene or CpG site
calculate_midpoint <- function(start, end) {
  return((start + end) / 2)
}

# Calculate midpoints for CAPN2 and each specified gene or CpG site
capn2_midpoint <- calculate_midpoint(ensembl_genes$start[ensembl_genes$external_gene_name == "CAPN2"],
                                     ensembl_genes$end[ensembl_genes$external_gene_name == "CAPN2"])

targets <- c("SCGB3A1", "GPR27", "ZNF542P", "GAL3ST2")

# Convert gene names to lowercase for comparison
ensembl_CpGsites$external_CpG_name <- tolower(ensembl_CpGsites$external_CpG_name)
targets <- tolower(targets)

distances_capn2 <- data.frame(main_gene = rep("CAPN2", length(targets)),
                        target = targets,
                        distance = NA)

for (i in seq_along(targets)) {
  target_start <- ensembl_CpGsites$start[ensembl_CpGsites$external_CpG_name == targets[i]]
  target_end <- ensembl_CpGsites$end[ensembl_CpGsites$external_CpG_name == targets[i]]
  
  # Check if the target exists in the CpG sites data frame
  if (length(target_start) > 0 && length(target_end) > 0) {
    target_midpoint <- calculate_midpoint(target_start, target_end)
    distances_capn2$distance[i] <- abs(capn2_midpoint - target_midpoint)
  } else {
    # Handle the case where the target is not found
    distances_capn2$distance[i] <- NA
  }
}

# Display the resulting data frame
print(distances_capn2)
# main_gene  target  distance
# 1     CAPN2 scgb3a1  43109429
# 2     CAPN2   gpr27 151947738
# 3     CAPN2 znf542p 167328343
# 4     CAPN2 gal3st2 223690430

write.table(distances_capn2, "distances_capn2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#####################################################
###Calculate distances between LRP1 and CpG sites###
#####################################################
# Convert start and end positions to numeric values
ensembl_genes$start <- as.numeric(ensembl_genes$start)
ensembl_genes$end <- as.numeric(ensembl_genes$end)

ensembl_CpGsites$start <- as.numeric(ensembl_CpGsites$start)
ensembl_CpGsites$end <- as.numeric(ensembl_CpGsites$end)

# Function to calculate the midpoint for a given gene or CpG site
calculate_midpoint <- function(start, end) {
  return((start + end) / 2)
}

# Calculate midpoints for LRP1 and each specified gene or CpG site
lrp1_midpoint <- calculate_midpoint(ensembl_genes$start[ensembl_genes$external_gene_name == "LRP1"],
                                     ensembl_genes$end[ensembl_genes$external_gene_name == "LRP1"])

targets <- c("BCOR", "TRIO")

# Convert gene names to lowercase for comparison
ensembl_CpGsites$external_CpG_name <- tolower(ensembl_CpGsites$external_CpG_name)
targets <- tolower(targets)

distances_lrp1 <- data.frame(main_gene = rep("LRP1", length(targets)),
                              target = targets,
                              distance = NA)

for (i in seq_along(targets)) {
  target_start <- ensembl_CpGsites$start[ensembl_CpGsites$external_CpG_name == targets[i]]
  target_end <- ensembl_CpGsites$end[ensembl_CpGsites$external_CpG_name == targets[i]]
  
  # Check if the target exists in the CpG sites data frame
  if (length(target_start) > 0 && length(target_end) > 0) {
    target_midpoint <- calculate_midpoint(target_start, target_end)
    distances_lrp1$distance[i] <- abs(lrp1_midpoint - target_midpoint)
  } else {
    # Handle the case where the target is not found
    distances_lrp1$distance[i] <- NA
  }
}

# Display the resulting data frame
print(distances_lrp1)
# main_gene target distance
# 1      LRP1   bcor 17059372
# 2      LRP1   trio 42985151

write.table(distances_lrp1, "distances_lrp1.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#####################################################
###Calculate distances between RPS6 and CpG sites###
#####################################################
# Convert start and end positions to numeric values
ensembl_genes$start <- as.numeric(ensembl_genes$start)
ensembl_genes$end <- as.numeric(ensembl_genes$end)

ensembl_CpGsites$start <- as.numeric(ensembl_CpGsites$start)
ensembl_CpGsites$end <- as.numeric(ensembl_CpGsites$end)

# Function to calculate the midpoint for a given gene or CpG site
calculate_midpoint <- function(start, end) {
  return((start + end) / 2)
}

# Calculate midpoints for RPS6 and each specified gene or CpG site
rps6_midpoint <- calculate_midpoint(ensembl_genes$start[ensembl_genes$external_gene_name == "RPS6"],
                                     ensembl_genes$end[ensembl_genes$external_gene_name == "RPS6"])

targets <- c("SRCIN1")

# Convert gene names to lowercase for comparison
ensembl_CpGsites$external_CpG_name <- tolower(ensembl_CpGsites$external_CpG_name)
targets <- tolower(targets)

distances_rps6 <- data.frame(main_gene = rep("RPS6", length(targets)),
                             target = targets,
                             distance = NA)

for (i in seq_along(targets)) {
  target_start <- ensembl_CpGsites$start[ensembl_CpGsites$external_CpG_name == targets[i]]
  target_end <- ensembl_CpGsites$end[ensembl_CpGsites$external_CpG_name == targets[i]]
  
  # Check if the target exists in the CpG sites data frame
  if (length(target_start) > 0 && length(target_end) > 0) {
    target_midpoint <- calculate_midpoint(target_start, target_end)
    distances_rps6$distance[i] <- abs(rps6_midpoint - target_midpoint)
  } else {
    # Handle the case where the target is not found
    distances_rps6$distance[i] <- NA
  }
}

# Display the resulting data frame
print(distances_rps6)
# main_gene target distance
# 1      RPS6 srcin1 16738987

write.table(distances_rps6, "distances_rps6.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#####################################################
###Calculate distances between RPS18 and CpG sites###
#####################################################
# Convert start and end positions to numeric values
ensembl_genes$start <- as.numeric(ensembl_genes$start)
ensembl_genes$end <- as.numeric(ensembl_genes$end)

ensembl_CpGsites$start <- as.numeric(ensembl_CpGsites$start)
ensembl_CpGsites$end <- as.numeric(ensembl_CpGsites$end)

# Function to calculate the midpoint for a given gene or CpG site
calculate_midpoint <- function(start, end) {
  return((start + end) / 2)
}

# Calculate midpoints for RPS18 and each specified gene or CpG site
rps18_midpoint <- calculate_midpoint(ensembl_genes$start[ensembl_genes$external_gene_name == "RPS18"],
                                     ensembl_genes$end[ensembl_genes$external_gene_name == "RPS18"])

targets <- c("TCEA3", "MYADML", "GMDS")

# Convert gene names to lowercase for comparison
ensembl_CpGsites$external_CpG_name <- tolower(ensembl_CpGsites$external_CpG_name)
targets <- tolower(targets)

distances_rps18 <- data.frame(main_gene = rep("RPS18", length(targets)),
                             target = targets,
                             distance = NA)

for (i in seq_along(targets)) {
  target_start <- ensembl_CpGsites$start[ensembl_CpGsites$external_CpG_name == targets[i]]
  target_end <- ensembl_CpGsites$end[ensembl_CpGsites$external_CpG_name == targets[i]]
  
  # Check if the target exists in the CpG sites data frame
  if (length(target_start) > 0 && length(target_end) > 0) {
    target_midpoint <- calculate_midpoint(target_start, target_end)
    distances_rps18$distance[i] <- abs(rps18_midpoint - target_midpoint)
  } else {
    # Handle the case where the target is not found
    distances_rps18$distance[i] <- NA
  }
}

# Display the resulting data frame
print(distances_rps18)
# main_gene target distance
# 1     RPS18  tcea3  4319796
# 2     RPS18 myadml 29213041
# 3     RPS18   gmds  2889783

write.table(distances_rps18, "distances_rps18.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#####################################################
###Calculate distances between EIF4G1 and CpG sites###
#####################################################
# Convert start and end positions to numeric values
ensembl_genes$start <- as.numeric(ensembl_genes$start)
ensembl_genes$end <- as.numeric(ensembl_genes$end)

ensembl_CpGsites$start <- as.numeric(ensembl_CpGsites$start)
ensembl_CpGsites$end <- as.numeric(ensembl_CpGsites$end)

# Function to calculate the midpoint for a given gene or CpG site
calculate_midpoint <- function(start, end) {
  return((start + end) / 2)
}

# Calculate midpoints for EIF4G1 and each specified gene or CpG site
eif4g1_midpoint <- calculate_midpoint(ensembl_genes$start[ensembl_genes$external_gene_name == "EIF4G1"],
                                     ensembl_genes$end[ensembl_genes$external_gene_name == "EIF4G1"])

targets <- c("PNMA8A", "COL6A6", "ZFHX4", "PLPPR3")

# Convert gene names to lowercase for comparison
ensembl_CpGsites$external_CpG_name <- tolower(ensembl_CpGsites$external_CpG_name)
targets <- tolower(targets)

distances_eif4g1 <- data.frame(main_gene = rep("EIF4G1", length(targets)),
                              target = targets,
                              distance = NA)

for (i in seq_along(targets)) {
  target_start <- ensembl_CpGsites$start[ensembl_CpGsites$external_CpG_name == targets[i]]
  target_end <- ensembl_CpGsites$end[ensembl_CpGsites$external_CpG_name == targets[i]]
  
  # Check if the target exists in the CpG sites data frame
  if (length(target_start) > 0 && length(target_end) > 0) {
    target_midpoint <- calculate_midpoint(target_start, target_end)
    distances_eif4g1$distance[i] <- abs(eif4g1_midpoint - target_midpoint)
  } else {
    # Handle the case where the target is not found
    distances_eif4g1$distance[i] <- NA
  }
}

# Display the resulting data frame
print(distances_eif4g1)
# main_gene target  distance
# 1    EIF4G1 pnma8a 137842087
# 2    EIF4G1 col6a6  53797318
# 3    EIF4G1  zfhx4 107633256
# 4    EIF4G1 plppr3 184312699

write.table(distances_eif4g1, "distances_eif4g1.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#####################################################
###Calculate distances between ACTB and CpG sites###
#####################################################
# Convert start and end positions to numeric values
ensembl_genes$start <- as.numeric(ensembl_genes$start)
ensembl_genes$end <- as.numeric(ensembl_genes$end)

ensembl_CpGsites$start <- as.numeric(ensembl_CpGsites$start)
ensembl_CpGsites$end <- as.numeric(ensembl_CpGsites$end)

# Function to calculate the midpoint for a given gene or CpG site
calculate_midpoint <- function(start, end) {
  return((start + end) / 2)
}

# Calculate midpoints for ACTB and each specified gene or CpG site
actb_midpoint <- calculate_midpoint(ensembl_genes$start[ensembl_genes$external_gene_name == "ACTB"],
                                     ensembl_genes$end[ensembl_genes$external_gene_name == "ACTB"])

targets <- c("PLPPR3")

# Convert gene names to lowercase for comparison
ensembl_CpGsites$external_CpG_name <- tolower(ensembl_CpGsites$external_CpG_name)
targets <- tolower(targets)

distances_actb <- data.frame(main_gene = rep("ACTB", length(targets)),
                               target = targets,
                               distance = NA)

for (i in seq_along(targets)) {
  target_start <- ensembl_CpGsites$start[ensembl_CpGsites$external_CpG_name == targets[i]]
  target_end <- ensembl_CpGsites$end[ensembl_CpGsites$external_CpG_name == targets[i]]
  
  # Check if the target exists in the CpG sites data frame
  if (length(target_start) > 0 && length(target_end) > 0) {
    target_midpoint <- calculate_midpoint(target_start, target_end)
    distances_actb$distance[i] <- abs(actb_midpoint - target_midpoint)
  } else {
    # Handle the case where the target is not found
    distances_actb$distance[i] <- NA
  }
}

# Display the resulting data frame
print(distances_actb)
# main_gene target distance
# 1      ACTB plppr3  5529806

write.table(distances_actb, "distances_actb.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#####################################################
###Calculate distances between DST and CpG sites###
#####################################################
# Convert start and end positions to numeric values
ensembl_genes$start <- as.numeric(ensembl_genes$start)
ensembl_genes$end <- as.numeric(ensembl_genes$end)

ensembl_CpGsites$start <- as.numeric(ensembl_CpGsites$start)
ensembl_CpGsites$end <- as.numeric(ensembl_CpGsites$end)

# Function to calculate the midpoint for a given gene or CpG site
calculate_midpoint <- function(start, end) {
  return((start + end) / 2)
}

# Calculate midpoints for DST and each specified gene or CpG site
dst_midpoint <- calculate_midpoint(ensembl_genes$start[ensembl_genes$external_gene_name == "DST"],
                                     ensembl_genes$end[ensembl_genes$external_gene_name == "DST"])

targets <- c("TRIO", "SFMBT1")

# Convert gene names to lowercase for comparison
ensembl_CpGsites$external_CpG_name <- tolower(ensembl_CpGsites$external_CpG_name)
targets <- tolower(targets)

distances_dst <- data.frame(main_gene = rep("DST", length(targets)),
                             target = targets,
                             distance = NA)

for (i in seq_along(targets)) {
  target_start <- ensembl_CpGsites$start[ensembl_CpGsites$external_CpG_name == targets[i]]
  target_end <- ensembl_CpGsites$end[ensembl_CpGsites$external_CpG_name == targets[i]]
  
  # Check if the target exists in the CpG sites data frame
  if (length(target_start) > 0 && length(target_end) > 0) {
    target_midpoint <- calculate_midpoint(target_start, target_end)
    distances_dst$distance[i] <- abs(dst_midpoint - target_midpoint)
  } else {
    # Handle the case where the target is not found
    distances_dst$distance[i] <- NA
  }
}

# Display the resulting data frame
print(distances_dst)
# main_gene target distance
# 1       DST   trio 42500655
# 2       DST sfmbt1  3596923

write.table(distances_dst, "distances_dst.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


#################################
###Merge all dataframes in one###
#################################

distances_genes_CpGs <- bind_rows(distances_actb, distances_capn2, distances_dst, distances_eif4g1, 
                                  distances_lrp1, distances_rps18, distances_rps6)

colnames(distances_genes_CpGs)[colnames(distances_genes_CpGs) == "distance"] <- "distance (bp)"

write.table(distances_genes_CpGs, "distance_between_genes_&_CpGs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
