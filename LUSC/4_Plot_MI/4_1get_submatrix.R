#!/usr/bin/env Rscript
library(tidyverse)
setwd("/home/mdiaz/lungsquamouscells/functional_enrichment")
########################PARAMETERS & PACKAGES
argms=commandArgs(trailingOnly=TRUE)
fun=argms[1]#"GO:"
subty="LUSC"

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
setwd("/home/mdiaz/lungsquamouscells/prepo_data/")
data=data.table::fread(paste(subty,"eigenNormi",sep='.'))
data=data[data$V1%in%features,]
setwd("~/lungsquamouscells/MI/submatrix/")
write_tsv(data,paste(fun,subty,"mtrx",sep='.'),col_names=F)#needed for puma
