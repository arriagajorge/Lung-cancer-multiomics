setwd("/home/mdiaz/lungsquamouscells/prepo_data/")

library(TCGAbiolinks)#2.20.1
library(biomaRt)#2.48.3  
library(tidyverse)

subtype=read.table("subtypeLUSC.tsv",header=T,sep='\t')
#get the data
mirnas <- GDCquery(project = "TCGA-LUSC",
                   data.category = "Transcriptome Profiling",
                   data.type = "miRNA Expression Quantification",
                   barcode=subtype$samples)
#Genome of reference: hg38
#https://api.gdc.cancer.gov/data/683367de-81c9-408c-85fd-2391f3e537ee
#says miRBase v.21 was used for harmonization annotation
GDCdownload(mirnas)
mir=GDCprepare(mirnas)
mirID=mir$miRNA_ID
write_tsv(as.data.frame(mirID),"miR.ids.map.tsv")

#####################################
#ESTE SCRIPT NO FUE UTILIZADO#########
#LOS miR.ids PUEDEN OBTENERSE EN:####
#https://mirbase.org/ftp/CURRENT/miRNA.xls.gz
