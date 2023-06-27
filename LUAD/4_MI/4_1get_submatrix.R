#!/usr/bin/env Rscript
library(tidyverse)
setwd("/home/jvasquez/Documents/Omics/LUAD/3_FuncEnrichment/")
########################PARAMETERS & PACKAGES
argms=commandArgs(trailingOnly=TRUE)
fun=argms[1]#"GO:0002709"
subty=argms[2]#"prox-inflam" #change names to files from prox.-inflam.... to prox-inflam.... 

#get the components linked to the function
if(length(grep("GO",fun))>0){
  enrich=read_tsv("BP.enrichment")
}else{
  enrich=read_tsv("KEGG.enrichment")	
}
comp=enrich%>%filter(ID==fun&subtype==subty)%>%
  dplyr::select(component)%>%unlist

#get the features selected in those components
selected=read_tsv(paste(subty,"stable",sep='.'))
features=selected%>%filter(component==unlist(comp))%>%
  distinct(variable)%>%unlist
print(paste("Function has",length(features),"features associated",sep=' '))

#get the data
data=data.table::fread(paste(subty,"eigenNormi",sep='.'))
data=data[data$V1%in%features,]
write_tsv(data,paste(fun,subty,"mtrx",sep='.'),col_names=F)#needed for puma