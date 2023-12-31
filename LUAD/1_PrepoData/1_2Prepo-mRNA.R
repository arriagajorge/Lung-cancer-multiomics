setwd("/home/jvasquez/Documents/Omics/LUAD/1_PrepoData/")
#FROM:
##  HANDS-ON:   NGS ANALYSIS  -->  Quantification files
## By Carlos Martínez Mira, Nov-2016
## By Sonia Tarazona, Nov-2016
##https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
## &&
## Data preparation: 
##      -Quality Control & bias removal
## By Cristobal Fresno - cristobalfresno@gmail.com
library(SummarizedExperiment)
library(TCGAbiolinks)
library(biomaRt)  
subtypeLUAD=read.table("subtypeLUAD.tsv",header=T, sep="\t")

xprssnLUAD <- GDCquery(project = "TCGA-LUAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts",
                       barcode=subtypeLUAD$samples)

df <- xprssnLUAD[[1]][[1]]
GDCdownload(xprssnLUAD)
expreLUADT <- GDCprepare(xprssnLUAD, summarizedExperiment = T)
expreLUAD <- GDCprepare(xprssnLUAD, summarizedExperiment = F)

tempLUAD = as.matrix(expreLUAD[,2:ncol(expreLUAD)])
rownames(tempLUAD) = expreLUAD$gene_id
expreLUAD = tempLUAD

stranded_firstLUAD <- expreLUADT@assays@data@listData[["stranded_first"]]
unstrandedLUAD <- expreLUADT@assays@data@listData[["unstranded"]]
stranded_secondLUAD <- expreLUADT@assays@data@listData[["stranded_second"]]
tpm_unstrandLUAD <- expreLUADT@assays@data@listData[["tpm_unstrand"]]
fpkm_unstrandLUAD <- expreLUADT@assays@data@listData[["fpkm_unstrand"]]
fpkm_uq_unstrandLUAD <- expreLUADT@assays@data@listData[["fpkm_uq_unstrand"]]
# "unstranded" "stranded_f" "stranded_s" "tpm_unstra" "fpkm_unstr" "fpkm_uq_un"

variables_ <- list(stranded_firstLUAD, unstrandedLUAD, stranded_secondLUAD, tpm_unstrandLUAD, 
                   fpkm_unstrandLUAD, fpkm_uq_unstrandLUAD)

colnames(stranded_firstLUAD) <- expreLUADT@colData@rownames
colnames(unstrandedLUAD) <- expreLUADT@colData@rownames
colnames(stranded_secondLUAD) <- expreLUADT@colData@rownames
colnames(tpm_unstrandLUAD) <- expreLUADT@colData@rownames
colnames(fpkm_unstrandLUAD) <- expreLUADT@colData@rownames
colnames(fpkm_uq_unstrandLUAD) <- expreLUADT@colData@rownames

write.table(fpkm_unstrandLUAD, "RNAseqLUAD2.tsv", sep = '\t', quote = F)
write.table(unstrandedLUAD, "RNAseqLUAD.tsv", sep = '\t', quote = F)
#stranded_firstLUAD <- cbind(gene_id = expreLUAD$gene_id[1:60660], stranded_firstLUAD)
#stranded_firstLUAD <- cbind(gene_name = expreLUAD$gene_name[1:60660], stranded_firstLUAD)
#stranded_firstLUAD <- cbind(gene_type = expreLUAD$gene_type[1:60660], stranded_firstLUAD)

# subtype to duplicates #one only
# i = substr(colnames(fpkm_unstrandLUAD), 1, 19)
i = substr(colnames(unstrandedLUAD), 1, 19)
j = i[duplicated(i)]
designExpLUAD=subtypeLUAD[c(which(!subtypeLUAD$samples%in%j),
                            as.numeric(sapply(which(subtypeLUAD$samples%in%j),rep,2))),]
designExpLUAD=designExpLUAD[order(match(designExpLUAD$samples,substr(colnames(expreLUAD),1,19))),]
designExpLUAD$barcode=colnames(fpkm_unstrandLUAD) # esto no cambia independientemente de con que transcriptoma estemos ocupando
write.table(designExpLUAD, "RNAdesignexp.tsv", sep = '\t', quote = F)
# expreLUAD[,"gene_type"][1:5]
# colnames(expreLUAD)[1:5]

add_gene_info <- function(dfLUAD){
  dfLUAD <- cbind(gene_type = expreLUAD[,"gene_type"][1:60660], dfLUAD)
  dfLUAD <- cbind(gene_name = expreLUAD[,"gene_name"][1:60660], dfLUAD)
  dfLUAD <- cbind(gene_id = rownames(expreLUAD)[1:60660], dfLUAD)
  return(dfLUAD)
}
unstrandedLUAD <- add_gene_info(unstrandedLUAD)
stranded_firstLUAD <- add_gene_info(stranded_firstLUAD)
stranded_secondLUAD <- add_gene_info(stranded_secondLUAD)
tpm_unstrandLUAD <- add_gene_info(tpm_unstrandLUAD)
fpkm_unstrandLUAD <- add_gene_info(fpkm_unstrandLUAD)
fpkm_uq_unstrandLUAD <- add_gene_info(fpkm_uq_unstrandLUAD)

dim(designExpLUAD)
# head(i); head(j); length(i); length(j); unique(j)

# keep only tenscript id not version numbers
rownames(unstrandedLUAD) <- unstrandedLUAD[,"gene_id"]
rownames(unstrandedLUAD) <- sapply(strsplit(rownames(unstrandedLUAD), ".", fixed=T),
                                   function(x) x[1])

rownames(stranded_firstLUAD) <- rownames(unstrandedLUAD)
rownames(stranded_secondLUAD) <- rownames(unstrandedLUAD)
rownames(tpm_unstrandLUAD) <- rownames(unstrandedLUAD)
rownames(fpkm_unstrandLUAD) <- rownames(unstrandedLUAD)
rownames(fpkm_uq_unstrandLUAD) <- rownames(unstrandedLUAD)

## "unstranded" "stranded_f" "stranded_s" "tpm_unstra" "fpkm_unstr" "fpkm_uq_un"
#annnotate GC content, length & biotype per transcript
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", 
                             "percentage_gene_gc_content", "gene_biotype",
                             "start_position","end_position","hgnc_id","hgnc_symbol"),
              filters = "ensembl_gene_id", 
              values=rownames(fpkm_unstrandLUAD),mart=mart) #its valid for every variable
myannot$length=abs(myannot$end_position-myannot$start_position)

#filter transcripts withouth annotation
myannot=myannot[myannot$gene_biotype=="protein_coding"&
                  myannot$hgnc_symbol!="",]
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
exprots_hgnc=unstrandedLUAD[rownames(fpkm_unstrandLUAD)%in%myannot$ensembl_gene_id,]
exprots_hgnc <- exprots_hgnc[,4:ncol(exprots_hgnc)]
dim(exprots_hgnc)
#exprots_hgnc[,"gene_id"]
#[1] 19400   201

##check duplicated probes
#myannot[myannot$hgnc_id == myannot$hgnc_id[duplicated(myannot$hgnc_id)],]
myannot2 <- myannot[unique(rownames(myannot)),]
dim(myannot2); dim(myannot) #same dimension
# > myannot$hgnc_id[duplicated(myannot$hgnc_id)]
# [1] "HGNC:30046" "HGNC:11582" "HGNC:33853" "HGNC:4876" 
which(myannot2$hgnc_id == "HGNC:30046"); which(myannot2$hgnc_id == "HGNC:11582")
which(myannot2$hgnc_id == "HGNC:33853"); which(myannot2$hgnc_id == "HGNC:4876")

# > which(!myannot$ensembl_gene_id%in%myannot3$ensembl_gene_id)
# [1] 10767
# ensembl_gene_id percentage_gene_gc_content   gene_biotype start_position end_position
# 41617 ENSG00000254093                      43.17 protein_coding       10764961     10839884
# hgnc_id hgnc_symbol length
# 41617 HGNC:30046       PINX1  74923

myannot2[c(18566,18728),]; myannot2[c(19327,19336),]
myannot2[c(12080,19340),]; myannot2[c(7581,19375),]

myannot3 <- myannot2[-c(18728,19327,19340,7581),]
dim(myannot2); dim(myannot3) #
# length(unique(rownames(stranded_firstLUAD)))
# length(unique(rownames(stranded_firstLUAD)))
# [1] 60616

#fpkm_LUAD <- fpkm_unstrandLUAD[unique(rownames(stranded_firstLUAD)),]

##################CHECK BIASES########################################################
library(NOISeq)
library(edgeR)

exprots_hgnc2 <- as.data.frame(exprots_hgnc[unique(rownames(exprots_hgnc)),],)
exprots_hgnc3 <- sapply(exprots_hgnc2, as.numeric)
rownames(exprots_hgnc3) <- rownames(exprots_hgnc[unique(rownames(exprots_hgnc)),])

#format data for noiseq
noiseqData = NOISeq::readData(data = exprots_hgnc3,
                              gc = myannot[,1:2],
                              biotype = myannot[,c(1,3)],factor=designExpLUAD,
                              length=myannot[,c(1,8)])
noiseqData2 = NOISeq::readData(data = exprots_hgnc3,
                               gc = myannot3[,1:2],
                               biotype = myannot3[,c(1,3)],factor=designExpLUAD,
                               length=myannot3[,c(1,8)])

#1)check expression bias per subtype
mycountsbio = NOISeq::dat(noiseqData, type = "countsbio", factor = "subtype")
# [1] "Warning: 249 features with 0 counts in all samples are to be removed for this analysis."
# [1] "Counts per million distributions are to be computed for:"
# [1] "normal"        "prox.-inflam"  "prox.-prolif." "TRU"
mycountsbio2 = NOISeq::dat(noiseqData2, type = "countsbio", factor = "subtype")


#patients with repeated measures
png("CountsOri.png")
explo.plot(mycountsbio2, plottype = "boxplot",samples = 1:1)
dev.off()
#2)check for low count genes
png("lowcountsOri.png")
explo.plot(mycountsbio2, plottype = "barplot", samples = 1:1)
dev.off()
png("lowCountThres.png")
hist(rowMeans(cpm(exprots_hgnc3,log=T)),ylab="genes",
     xlab="mean of log CPM",col="gray")
abline(v=0,col="red")
dev.off()

#3)check for transcript composition bias
#each sample s is compared to a reference r (which can be arbitrarily chosen).
#by computing M values=log2(countss = countsr). 
#Confidence intervals for the M median is computed by bootstrapping.
#If the median of M values for each comparison is not in the CI, the deviation
# of the sample is significant, therefore, normalization is needed 
mycd = NOISeq::dat(noiseqData2, type = "cd", norm = FALSE) #slooooow
# [1] "Warning: 250 features with 0 counts in all samples are to be removed for this analysis."
# [1] "Reference sample is: TCGA-64-5781-01A-01R-1628-07"
# [1] "Diagnostic test: FAILED. Normalization is required to correct this bias."
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   165     40
png("MvaluesOri.png")
explo.plot(mycd,samples=sample(1:ncol(exprots_hgnc3),10))
dev.off()

#4)check for length & GC bias
#A cubic spline regression model is fitted. Both the model p-value and the coefficient
# of determination (R2) are shown. If the model p-value is significant and R2 value is
# high (more than 70%) the exp,ression depends on the feature
myGCcontent <- NOISeq::dat(noiseqData2, type = "GCbias", factor = "subtype")
png("GCbiasOri.png",width=1000)
par(mfrow=c(1,1))
sapply(1:1, function(x) explo.plot(myGCcontent, samples = x))
dev.off()
#The GC-content of each gene does not change from sample to sample, so it can be expected to
#have little effect on differential expression analyses to a first approximation
mylenBias <- NOISeq::dat(noiseqData2, k = 0, type = "lengthbias",
                         factor = "subtype")
png("lengthbiasOri.png",width=1000)
par(mfrow=c(1,1))
sapply(1:1,function(x) explo.plot(mylenBias, samples = x))
dev.off()
#BUT, since the gene has the same length in all your samples, there is no need to divide by the gene length

# #5) check for batch effect
# myPCA = NOISeq::dat(noiseqData2, type = "PCA", norm = F, logtransf = F)
# png("PCA_Ori.png")
# explo.plot(myPCA, samples = c(1,2), plottype = "scores",
#            factor = "subtype")
# dev.off()

#################SOLVE BIASES######################################################
library(EDASeq)

#1) filter low count genes.
#CPM=(counts/fragments sequenced)*one million.
#Filtering those genes with average CPM below 1, would be different
#to filtering by those with average counts below 1.
countMatrixFiltered = filtered.data(exprots_hgnc3, factor = "subtype",
                                    norm = FALSE, depth = NULL, method = 1, cpm = 0, p.adj = "fdr")
#10943 features are to be kept for differential expression analysis with filtering method 1
myannot3=myannot3[myannot3$ensembl_gene_id%in%rownames(countMatrixFiltered),]
myannot=myannot[myannot$ensembl_gene_id%in%rownames(countMatrixFiltered),]


#all names must match
mydataEDA <- newSeqExpressionSet(
  counts=as.matrix(countMatrixFiltered),
  featureData=data.frame(myannot3,row.names=myannot3$ensembl_gene_id),
  phenoData=data.frame(designExpLUAD,row.names=designExpLUAD$barcode))

mydataEDA2 <- newSeqExpressionSet(
  counts=as.matrix(countMatrixFiltered),
  featureData=data.frame(myannot,row.names=myannot$ensembl_gene_id),
  phenoData=data.frame(designExpLUAD,row.names=designExpLUAD$barcode))
#order for less bias
gcFull <- withinLaneNormalization(mydataEDA2, 
                                  "percentage_gene_gc_content", which = "full")#corrects GC bias 
lFull <- withinLaneNormalization(gcFull, "length", which = "full")#corrects length bias 
fullfullTMM <- NOISeq::tmm(normCounts(lFull), long = 1000, lc = 0, k = 0)
#norm.counts <- betweenLaneNormalization(normCounts(lFull),
# which = "median", offset = FALSE)
#FAILED PASSED 
#   290    518
noiseqData = NOISeq::readData(data = fullfullTMM, factors=designExpLUAD)
#cd has to preceed ARSyN or won't work
mycd=NOISeq::dat(noiseqData,type="cd",norm=TRUE)
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"]) #sometimes change values
# FAILED PASSED
#     0     195
# "Diagnostic test: PASSED."
# we only have 1 subtype
#############################SOLVE BATCH EFFECT#######################################################
# myPCA = NOISeq::dat(noiseqData, type = "PCA", norm = T, logtransf = F)
# png("preArsyn.png")
# explo.plot(myPCA, samples = c(1,2), plottype = "scores",
#            factor = "subtype")
# dev.off()
# ffTMMARSyn=ARSyNseq(noiseqData, batch = F, factor = "subtype",
#                     norm = "n",  logtransf = T)
# myPCA = NOISeq::dat(ffTMMARSyn, type = "PCA", norm = T,logtransf = T)
# png("postArsyn.png")
# explo.plot(myPCA, samples = c(1,2), plottype = "scores", 
#            factor = "subtype")
# dev.off()
# 
# explo.plot(myPCA, samples = c(1,2), plottype = "scores",
#            factor = "race")
# explo.plot(myPCA, samples = c(1,2), plottype = "scores",
#            factor = "gender")
#nothing relevant in the last two pca
#############################FINAL QUALITY CHECK#######################################################
#change fullfullTMM
noiseqData <- NOISeq::readData(data = fullfullTMM, gc = myannot[,1:2],
                               biotype = myannot[,c(1,3)],factor=designExpLUAD,
                               length=myannot[,c(1,8)])
mycountsbio = NOISeq::dat(noiseqData, type = "countsbio", factor = "subtype",
                          norm=T)
png("CountsFinal.png")
explo.plot(mycountsbio, plottype = "boxplot",samples=1:1)
dev.off()
myGCcontent <- NOISeq::dat(noiseqData, k = 0, type = "GCbias", 
                           factor = "subtype",norm=T)
png("GCbiasFinal.png",width=1000)
par(mfrow=c(1,1))
sapply(1:1,function(x) explo.plot(myGCcontent, samples = x))
dev.off()
mylenBias <- NOISeq::dat(noiseqData, k = 0, type = "lengthbias", 
                         factor = "subtype",norm=T)
png("lengthbiasFinal.png",width=1000)
par(mfrow=c(1,1))
sapply(1:1,function(x) explo.plot(mylenBias, samples = x))
dev.off()

#############################RESOLVE DUPLICATES & SAVE##################################################
#get duplicates
i=designExpLUAD$samples[duplicated(designExpLUAD$samples)]
#get sample barcode per sample
i=lapply(i,function(x) designExpLUAD$barcode[designExpLUAD$samples==x])
#separate duplicates
final=fullfullTMM
duplis=final[,colnames(final)%in%unlist(i)]
prefi=final[,!colnames(final)%in%unlist(i)]
#average duplicates
temp=do.call(cbind,lapply(i,function(x) 
  rowMeans(duplis[,colnames(duplis)%in%x])))
#identify samples with barcode 
colnames(temp)=designExpLUAD$samples[duplicated(designExpLUAD$samples)]
colnames(prefi)=substr(colnames(prefi),1,19)
#joint matrices
final=cbind(prefi,temp)
dim(final)
# [1] 11399   193
final=final[,order(match(colnames(final),subtypeLUAD$samples))]

#in case you have zeros, in general if you use TPM, FPKM ...
# eliminarCerosbyRows <- function(df){
#   vectorCeros = c()
#   for (i in 1:dim(df)[1]) {
#     if (df[i,1]==0) {
#       if(sum(df[i,]) == 0){
#         vectorCeros <- c(vectorCeros, i)
#       }
#     }
#   }
#   return(vectorCeros)
# }
# final2 <- final[-eliminarCerosbyRows(final),]
# dim(final2)

write.table(final,"RNAseqnormalized.tsv",sep='\t',quote=F)
