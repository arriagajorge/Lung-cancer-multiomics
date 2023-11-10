setwd("/home/jvasquez/Documents/Lung-Cancer/LUAD/prepo")
library(tidyverse)
library(biomaRt)
#downloaded from https://support.illumina.com/downloads/humanmethylation450_15017482_v1-2_product_files.html
methy=read_csv("humanmethylation450_15017482_v1-2.csv",skip=7) 
#drop control probes
methy=methy[grep('^[0-9]',methy$IlmnID,perl=T,invert=T),]
#only keep probe ID, mapping genes & affected position
methy=methy%>%dplyr::select(IlmnID,UCSC_RefGene_Accession,
                            UCSC_RefGene_Group)%>%filter(!is.na(UCSC_RefGene_Accession))%>%
  separate_rows(c(UCSC_RefGene_Accession,UCSC_RefGene_Group),
                sep=';',convert=T)
#add entrez id, needed for most enrichment tools	
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("entrezgene_id",
                             "refseq_mrna","refseq_ncrna"), mart=mart)
myannot=myannot%>%pivot_longer(-1,names_to="type",values_to="refseq")

colnames(methy)[2] <- "refseq"
# methy=myannot%>%select(entrezgene_id,refseq)%>%
#   merge(methy,by="refseq",all.y=T)

myannot_subset <- subset(myannot, select = c("entrezgene_id", "refseq"))
merged_data <- merge(myannot_subset, methy, by = "refseq", all.y = TRUE)

methy = merged_data
length(unique(methy$IlmnID[is.na(methy$entrezgene_id)]))
#[1] 27106
#[1] 26912 out of 365860 probes don't map to entrez
write_tsv(methy,"MapMethy.tsv")
#https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/infinium-methylationepic-manifest-column-headings.pdf
