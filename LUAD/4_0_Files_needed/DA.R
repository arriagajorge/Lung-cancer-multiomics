setwd("/home/jvasquez/Documents/Omics/LUAD/1_PrepoData/")
library(limma)
library(data.table)
library(ggplot2)
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#https://f1000research.com/articles/5-1408/v1

#################RNA###########################################
subtype=read.table("subtypeLUAD.tsv",header=T,sep='\t')
expre=fread("RNAseqnormalized.tsv")
expre=as.matrix(expre[,2:ncol(expre)],rownames=expre$V1)

#set comparisons
#~0 gives a model where each coefficient corresponds to a group mean
subtype$subtype=factor(subtype$subtype)
design=model.matrix(~0+subtype,subtype)
colnames(design) <- c("LUAD", "normal")
#fix names
#colnames(design)=gsub("subtype","",colnames(design))
contr.mtrx=makeContrasts(
  LUAD_normal=LUAD-normal,
  levels=design)

#check if count&variance are indi
#if counts are more variable at lower expression, voom makes the 
#data “normal enough”, 
#wont work with ur log2-normalized stuff
#v=voom(expre,design,plot=T,save.plot=T)#no need if fitted curve is smooth 
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

#fit a linear model using weighted least squares for each gene
fit=lmFit(log2(expre),design)#lots of NA, if no log2 lfc is huge & u get lot of DE.genes
fitSubtype = contrasts.fit(fit, contr.mtrx)
#treat is better than fc+p.val thresholds, that increase FP
tfitSubtype=treat(fitSubtype, lfc = log2(1.5))
#log2(1.2 or 1.5) will usually give DE genes with fc => 2
#depending on the sample size and precision of the experiment
DE.genes=lapply(1:1,function(x) 
  topTreat(tfitSubtype,coef=x,n=nrow(expre)))
names(DE.genes)=colnames(contr.mtrx)
sapply(DE.genes,function(x) sum(x$adj.P.Val<0.05))
#        5334         4455         4075         4982 
##`pdf("DEgenes.pdf")`
#par(mfrow=c(2,2))
#sapply(1:4,function(x) plotMA(fitSubtype,coef=x))
#sapply(1:4,function(x) {
#	volcanoplot(tfitSubtype,coef=x)
#	abline(h=-log2(1.2),col="red")
#})
#dev.off()
temp=do.call(rbind,lapply(1:1,function(x) 
  cbind(contrast=names(DE.genes)[x],
        ensembl_gene_id=rownames(DE.genes[[x]]),
        DE.genes[[x]])))
png("logFC.png")
ggplot(temp,aes(y=logFC,x=contrast,color=contrast))+
  geom_boxplot()+theme(legend.position="none")
dev.off()
write.table(temp,"DE.genes.tsv",sep='\t',quote=F,row.names=F)
#next:GSEA
#################miRNA###########################################
mir=fread("miRNAseqNormi.tsv")
mir=as.matrix(mir[,2:ncol(mir)],rownames=mir$V1)

v=voom(mir,design,plot=T,save.plot=T)#coz mir is normalized, but no log2 transformed
fit=lmFit(v,design)
fitSubtype = contrasts.fit(fit, contr.mtrx)
#treat doesn't find any DE miR
#tfitSubtype=treat(fitSubtype, lfc = log2(1.2))
#DE.miR=lapply(1:4,function(x) 
#	topTreat(tfitSubtype,coef=x,n=nrow(mir)))
temp=eBayes(fitSubtype)
DE.miR=lapply(1:1,function(x) 
  topTable(temp,coef=x,n=Inf))
names(DE.miR)=colnames(contr.mtrx)
sapply(DE.miR,function(x) sum(x$adj.P.Val<0.05))
#Basal_Normal  
#         47
temp=do.call(rbind,lapply(1:1,function(x) 
  cbind(contrast=names(DE.miR)[x],
        id=rownames(DE.miR[[x]]),
        DE.miR[[x]])))
write.table(temp,"DE.miR.tsv",sep='\t',quote=F,row.names=F)

#################METHYLATION###########################################
#get covariates
clin <- TCGAbiolinks::GDCquery_clinic("TCGA-LUAD","clinical")
colnames(subtype)[3]="submitter_id"
subtype=merge(subtype,clin,by="submitter_id",all.x=T)
#apply(subtype,2,table) choose cols with various categories
subtype$ajc
subtype=subtype[,c("submitter_id",
                   "samples",
                   "tissue",
                   "subtype",
                   "age_at_diagnosis",
                   "treatments_pharmaceutical_treatment_or_therapy",
                   "treatments_radiation_treatment_or_therapy",
                   "race.x",
                   "ajcc_pathologic_m",
                   "ajcc_pathologic_n",
                   "morphology",
                   "ajcc_pathologic_t",
                   "prior_malignancy",
                   "primary_diagnosis",
                   "ajcc_pathologic_stage.x")]
write.table(subtype,"subtypeDA.tsv",sep='\t',quote=F,row.names=F)
methy=data.table::fread("methyM.tsv")
methy=as.matrix(methy[,2:ncol(methy)],rownames=methy$V1)
subtype=subtype[order(match(subtype$samples,colnames(methy))),]
#change category names for numbers or makeContrast wont work 
subtype=cbind(subtype[,1:5],
              apply(subtype[,6:15],2,function(x)
                factor(as.numeric(factor(x)))))
subtype[is.na(subtype)]=0#force NA to be another category or ur loose samples

#GIVEN THE LARGE N(>>5) AND CANCER EFFECT, ALL METHODS ARE EXPECTED TO WORK SIMILARLY

#1st approach: limma (control for covariates)
design=model.matrix(~0+subtype+age_at_diagnosis+morphology+
                      ajcc_pathologic_m+ajcc_pathologic_n+ajcc_pathologic_t+
                      primary_diagnosis+ajcc_pathologic_stage.x,subtype)
#used this complex model coz methy didn't went through ARSYN 
#no treatment covarites/prior_malignancy coz groups're mixed in PCA
colnames(design)=gsub("subtype","",colnames(design))
contr.mtrx=makeContrasts(
  LUAD_normal=LUAD-normal,
  levels=design)
# contr.mtrx=makeContrasts(
#   Basal_Normal=Basal-Normal,
#   Her2_Normal=Her2-Normal,
#   LumA_Normal=LumA-Normal,
#   LumB_Normal=LumB-Normal,
#   levels=design)
fit=lmFit(methy,design)
#Coefficients not estimable: primary_diagnosis
#Partial NA coefficients for 393132 probe(s) 
fitSubtype = contrasts.fit(fit, contr.mtrx)
#use treat to get lower dm cpgs than eBayes would
tfitSubtype=treat(fitSubtype, lfc = log2(1.5))
DMcpgs=lapply(1:1,function(x) 
  topTreat(tfitSubtype,coef=x,n=Inf))
names(DMcpgs)=colnames(contr.mtrx)
sapply(DMcpgs,function(x) sum(x$adj.P.Val<0.01))
# LUAD_normal 
# 1899

#when desing only considers subtype, DMs don't change much
#Basal_Normal  Her2_Normal  LumA_Normal  LumB_Normal 
#       64687        79396        88830       115412 
#but Maksimovic2015 says unadjusted limma has lots of FP
#is this approach enriched on FP??????
temp=do.call(rbind,lapply(1:1,function(x) 
  cbind(contrast=names(DMcpgs)[x],
        id=rownames(DMcpgs[[x]]),
        DMcpgs[[x]])))
write.table(temp,"DMcpgs.tsv",sep='\t',quote=F,row.names=F)

#2nd apporach; RUV-inverse (outperforms existing methods in DA of 450k data)
#stage 1: find empirical controls not associated with subtypes 
design=model.matrix(~0+subtype,subtype)
colnames(design)=gsub("subtype","",colnames(design))
contr.mtrx=makeContrasts(
  LUAD_normal = LUAD - normal,
  levels=design)
# 
# contr.mtrx=makeContrasts(
#   Basal_Normal=Basal-Normal,
#   Her2_Normal=Her2-Normal,
#   LumA_Normal=LumA-Normal,
#   LumB_Normal=LumB-Normal,
#   levels=design)
fit=lmFit(methy,design)
fitSubtype = contrasts.fit(fit, contr.mtrx)
temp=eBayes(fitSubtype) 
DMcpgs=lapply(1:1,function(x) topTable(temp,coef=x,n=Inf))

#Basal_Normal  Her2_Normal  LumA_Normal  LumB_Normal 
#      246810       228203       261208       256139 
ctlimma=lapply(DMcpgs,function(x) x$adj.P.Val > 0.5)
sapply(ctlimma,table)
# [,1]
# FALSE 194836
# TRUE  219831

#stage 2: differential methylation
library(missMethyl)
# rfit <- sapply(1:3,function(x)
#   #x=1
#   #compare each subtype with the normal tissue (4th unique subtype) 
#   RUVfit(Y=methy[,subtype$subtype%in%unique(subtype$subtype)[c(x,4)]],
#          X=factor(subtype$subtype[subtype$subtype%in%unique(subtype$subtype)[c(x,4)]]),
#          ctl=ctlimma[[x]]))


# fRUVfit <- function(x){
#   return(RUVfit(Y=methy[,subtype$subtype%in%unique(subtype$subtype)[c(x,2)]],
#                 X=factor(subtype$subtype[subtype$subtype%in%unique(subtype$subtype)[c(x,2)]]),
#                 ctl=ctlimma[[x]]))
# }
x=1
Y=methy[,subtype$subtype%in%unique(subtype$subtype)[c(x,2)]]
X=factor(subtype$subtype[subtype$subtype%in%unique(subtype$subtype)[c(x,2)]])
ctl=ctlimma[[x]]

rf1 <- RUVfit(Y=Y, X=X, ctl=ctl)
#; rf2 <- fRUVfit(2); rf3 <- fRUVfit(3); #fails with lapply

#adjust for unwanted variation #same fails with lapply
r1 <- RUVadj(Y=methy[,subtype$subtype%in%unique(subtype$subtype)[c(1,2)]], fit=rf1)
# 
# r2 <- RUVadj(Y=methy[,subtype$subtype%in%unique(subtype$subtype)[c(2,4)]], rf2)
# r3 <- RUVadj(Y=methy[,subtype$subtype%in%unique(subtype$subtype)[c(3,4)]], fit=rf3)
rfit <- list(r1)#, r2, r3)

#get the statistics of each cpg
DMcpg=BiocGenerics::lapply(rfit,function(x) topRUV(x,num=Inf))
names(DMcpgs)=colnames(contr.mtrx)
sapply(DMcpgs,function(x) sum(x$F.p.BH<0.05))
sapply(DMcpgs,function(x) sum(x$adj.P.Val<0.05))
# LUAD_normal 
# 44355 
temp=do.call(rbind,lapply(1:1,function(x) 
  cbind(contrast=names(DMcpgs)[x],
        id=rownames(DMcpgs[[x]]),
        DMcpgs[[x]])))
write.table(temp,"DMcpgs-RUV.tsv",sep='\t',quote=F,row.names=F)