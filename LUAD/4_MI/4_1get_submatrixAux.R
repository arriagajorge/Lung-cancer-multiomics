#!/usr/bin/env Rscript
library(tidyverse)
setwd("/home/jvasquez/Documents/Omics/LUAD/3_FuncEnrichment/")
bp <- read_tsv("BP.enrichment")
kegg <- read_tsv("KEGG.enrichment")

for (i in 1:dim(bp)[1]) {
  fun <- as.character(bp[i, "ID"])
  subty <- as.character(bp[i,"subtype"])
  
  command <-  paste("Rscript 4_1get_submatrix.R", fun, subty)
  setwd("/home/jvasquez/Documents/Omics/LUAD/4_MI/")
  system(command)
  setwd("/home/jvasquez/Documents/Omics/LUAD/3_FuncEnrichment/")
  arch <- paste(fun,subty,"mtrx",sep='.')
  system(paste(paste("mv ", arch), " /home/jvasquez/Documents/Aracne/ARACNE-multicore/launch"))
  print(i)
}

for (i in 1:dim(kegg)[1]) {
  fun <- as.character(kegg[i, "ID"])
  subty <- as.character(kegg[i,"subtype"])
  
  command <-  paste("Rscript 4_1get_submatrix.R", fun, subty)
  setwd("/home/jvasquez/Documents/Omics/LUAD/4_MI/")
  system(command)
  setwd("/home/jvasquez/Documents/Omics/LUAD/3_FuncEnrichment/")
  arch <- paste(fun,subty,"mtrx",sep='.')
  system(paste(paste("mv ", arch), " /home/jvasquez/Documents/Aracne/ARACNE-multicore/launch"))
  print(i)
}