setwd("/home/mdiaz/lungsquamouscells/sgcca/selected/")
library(tidyverse)

# subtype=commandArgs(trailingOnly=TRUE)
subtype="LUSC" 

files=list.files()
files=files[grep(subtype,files)]
subsa=lapply(files,read_tsv)
final=read_tsv(paste("../",subtype,".selected",sep=''))

subsa=lapply(subsa,function(x) x%>%dplyr::count(feature))
subsa=reduce(subsa,full_join,by="feature")
colnames(subsa)[2:101]=paste("r",1:100,sep="")
final=final%>%count(variable)
colnames(final)[1]="feature"
subsa=merge(final,subsa,by="feature",all.x=T)

write_tsv(subsa,paste(subtype,"subsampled",sep='.'))
#temp=apply(subsa[,2:101],1,table)