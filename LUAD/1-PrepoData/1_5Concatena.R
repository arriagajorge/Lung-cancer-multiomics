setwd("/home/jvasquez/Documents/Omics/LUAD/1_PrepoData/")
#!/usr/bin/env Rscript
library(data.table)
subtypeLUAD=read.table("subtypeLUAD.tsv",header=T,sep='\t')
# revert comments for unnormalized data
# expre=fread("RNAseqnormalized.tsv")
# miR=fread("miRNAseqNormi.tsv")
# methy=fread("methyM.tsv")
system("mkdir multiomics")
system("cp RNAseqnormalized.tsv ./multiomics")
system("cp miRNAseqNormi.tsv ./multiomics")
system("cp methyM.tsv ./multiomics")
files=list.files("./multiomics",full.names=T)
#files=files[grep("eigenNormi",files)]
data=lapply(files,fread)
#expre=as.matrix(expre[,2:ncol(expre)],rownames=expre$V1)
#miR=as.matrix(miR[,2:ncol(miR)],rownames=miR$V1)
#methy=as.matrix(methy[,2:ncol(methy)],rownames=methy$V1)

data=lapply(data,function(x) as.matrix(x[,2:ncol(x)],rownames=x$V1))
names(data)=gsub("./multiomics/","",files)
names(data)=gsub(".tsv","",names(data))
names(data)=gsub("RNAseqnormalized","transcripts",
                 gsub("miRNAseqNormi","miRNAs",gsub("methyM","CpGs",names(data))))
print(sapply(data,dim))
# CpGs miRNAs transcripts
# [1,] 414667    241       11475
# [2,]    188    188         188
#print(sapply(data,function(x) head(rownames(x))))

#choose methy order
#subtype=subtype[order(match(subtype$samples,colnames(methy))),]
#expre=expre[,order(match(colnames(expre),subtype$samples))]
#miR=miR[,order(match(colnames(miR),subtype$samples))]
subtype=subtypeLUAD[order(match(subtypeLUAD$samples,colnames(data$CpGs))),]
data[2:3]=lapply(data[2:3],function(x)
  x[,order(match(colnames(x),subtype$samples))])
print(names(data))
#"CpGs"        "miRNAseq"    "transcripts"

subtypeLUAD$subtype <- as.factor(subtypeLUAD$subtype)
#data per subtype
concatenated=lapply(levels(subtypeLUAD$subtype),function(x) 
  list(CpGs=data$CpGs[,subtypeLUAD$subtype==x],
       transcripts=data$transcripts[,subtypeLUAD$subtype==x],
       miRNA=data$miRNAs[,subtypeLUAD$subtype==x]))

names(concatenated)=levels(subtypeLUAD$subtype)


##########################################matrix per subtype
concatenated=lapply(concatenated,function(x) do.call(rbind,x))
print(sapply(concatenated,dim))
# LUAD
# [1,] 426383
# [2,]    188
lapply(1:1,function(x) write.table(concatenated[[x]],
                                   paste(names(concatenated)[x],"mtrx",sep='.'),sep='\t',quote=F))