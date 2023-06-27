setwd("/home/jvasquez/Documents/Omics/LUAD/1_PrepoData/")
library(TCGAbiolinks)
library(data.table)

subtypeLUAD=read.table("subtypeLUAD.tsv",header=T,sep='\t')
#get the data
# mthyltn2 <-  GDCquery(project = "TCGA-LUAD",
#                      data.category = "DNA Methylation",
#                      platform="Illumina Human Methylation 450",
#                      barcode=subtypeLUAD$samples)
# Error in GDCdownload(mthyltn) : 
#   We can only download one data type. Please use data.type argument in GDCquery to filter results.

mthyltn <-  GDCquery(project = "TCGA-LUAD",
                     data.category = "DNA Methylation",
                     data.type = "Methylation Beta Value",
                     platform="Illumina Human Methylation 450",
                     barcode=subtypeLUAD$samples)

# mthyltn <-  GDCquery(project = "TCGA-LUAD",
#                      data.category = "DNA Methylation",
#                      data.type = "Masked Intensities",
#                      platform="Illumina Human Methylation 450",
#                      barcode=subtypeLUAD$samples)
# In addition: Warning messages:
#1: In fread(x, select = 2, stringsAsFactors = F) :
#  Previous fread() session was not cleaned up properly. Cleaned up ok at the beginning of this fread() call.
#2: In fread(x, select = 2, stringsAsFactors = F) :

GDCdownload(mthyltn)
files=list.files("GDCdata",full.names=T,recursive=T)
#files=list.files("GDCdata",full.names=T,recursive=T)
#GDCdownload will download 197 files. A total of X GB

###############################################
#methy=do.call(cbind,pbapply::pbsapply(files,function(x) 
#  fread(x,select=2,stringsAsFactors=F)))

library(sesameData)
library(sesame)
mthyltn2=GDCprepare(mthyltn)
methy <- mthyltn2@assays@data[[1]]
dim(methy)
#[1] 485577    197
write.table(methy,"methy.tsv",sep='\t',quote=F)

#
# methy_2=BiocGenerics::do.call(cbind,pbapply::pbsapply(files,function(x) 
#   fread(x,select=2,stringsAsFactors=F)))
#######drop probes with too many na########################################
#before noise-prone filter or the process will be slooow

#subtype to duplicate
i=substr(colnames(methy),1,19)
j=i[duplicated(i)]
designMethy=subtypeLUAD[c(which(!subtypeLUAD$samples%in%j),
                          as.numeric(sapply(which(subtypeLUAD$samples%in%j),rep,2))),1:7]
#needed coz names are not equal to expression data but barcodes do
designMethy$barcode=unlist(sapply(designMethy$samples,function(x)
  colnames(methy)[i==x][1]))
designMethy$barcode[designMethy$samples%in%j]=rev(colnames(methy)[
  which(i%in%j)])
designMethy=designMethy[order(match(designMethy$barcode,
                                    colnames(methy))),]
total=table(subtypeLUAD$subtype)
#NA per subtype
nas=lapply(names(total),function(x) 
  rowSums(is.na(methy[,designMethy$subtype==x])))
#keep probes with NA in less than 25% of samples of all subtypes
i=unique(unlist(lapply(1:1,function(x) which(nas[[x]]<total[x]*.25))))
methy=methy[i,]
dim(methy)
#[1] 417031    844
#[1] 416370    398
############filter noise-prone probes##############################
#get probes annotation
annot=fread(files[1])
colnames(annot)[1]="probe"
annot=annot[annot$probe%in%rownames(methy),]
annot=annot[order(match(annot$probe,rownames(methy))),]
#sum(annot$probe==rownames(methy)) #order is the same

#?????
#sex chrs should be dropped if mixed gender
#since this is not the case, keep chrX
#https://bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
#methy=methy[!annot$Chromosome%in%c("chrX","chrY"),]
#dim(methy)

#drop probes with ambigous mapping
#https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Methylation_LO_Pipeline/
methy_=methy[annot$Chromosome!="*",]
nrow(methy_)
#[1] 417031

#polymorphisms may affect DNAm measurements
methy_grep=methy[grep("rs",rownames(methy),invert=T),]
nrow(methy_grep)
#[1] 417031 all values for elements that do not match.

#######impute missing data########################################
library(doParallel)
library(impute)

#separate per subtype
methy=lapply(names(total),function(x) 
  methy[,designMethy$subtype==x])
nas=lapply(methy,function(x) rowSums(is.na(x)))
sapply(nas,max)
#[1]  46
sapply(nas,function(x) sum(x>0))#cell to estimate
#[1] 42825

#impute missing values
no_cores=detectCores()-1
registerDoParallel(cores=no_cores)
#cl=makeCluster(no_cores,type="FORK")
#takes a looong time
#methy=parLapply(cl,methy,function(x) 
#  impute.knn(x,k=15,maxp=nrow(x))$data)#maxp forces impute, even if too many NAs
#Buuren and Groothuis suggested that 15 variables are generally 
#sufficient for imputation
#stopCluster(cl)
#alternative
methy=pbapply::pblapply(methy,function(x) impute.knn(x,k=15,maxp=nrow(x))$data)

write.table(do.call(cbind,methy),"methyNormi.tsv",sep='\t',quote=F)
methydf <- read.table("methyNormi.tsv",header=T,sep='\t')

#######beta to M values########################################
#Beta values have severe heteroscedasticity for highly methylated or
#unmethylated CpG sites. M-values provide much better performance in
#terms of detection rate and true positive rate for both highly 
#methylated and unmethylated CpG sites
library(ggplot2)

#check beta distributions
temp=lapply(methy,function(x) sample(x,10000))
temp=as.data.frame(do.call(rbind,lapply(1:1,function(x) 
  cbind(temp[[x]],levels(as.factor(designMethy$subtype))[x]))))
colnames(temp)=c("beta","subtype")
temp$beta=as.numeric(as.character(temp$beta))
png("betaDistr.png")
ggplot(temp,aes(x=beta))+
  geom_density(aes(fill=subtype,color=subtype,y=..scaled..),
               alpha=0.15)+theme(text=element_text(size=18))
dev.off()

#transform beta to M as in:
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
mval=function(beta){log2(beta/(1-beta))}
m=pbapply::pblapply(methy,function(x) apply(x,c(1,2),mval))
#check M distributions
temp=lapply(m,function(x) sample(x,10000))
temp=as.data.frame(do.call(rbind,lapply(1:1,function(x) 
  cbind(temp[[x]],levels(as.factor(designMethy$subtype))[x]))))
colnames(temp)=c("M","subtype")
temp$M=as.numeric(as.character(temp$M))
png("MDistr.png")
ggplot(temp,aes(x=M))+
  geom_density(aes(fill=subtype,color=subtype,y=..scaled..),
               alpha=0.1)+theme(text=element_text(size=18))
dev.off()
# Warning message: The dot-dot notation (`..scaled..`) was deprecated in ggplot2 3.4.0.
# ℹ Please use `after_stat(scaled)` instead.
# This warning is displayed once every 8 hours.
# Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated. 
#######average duplicates########################################
m=do.call(cbind,m)
i=substr(colnames(m),1,19)
j=i[duplicated(i)]
prefi=m[,!i%in%j]
duplis=m[,i%in%j]

temp=do.call(cbind,lapply(j,function(x) 
  rowMeans(duplis[,substr(colnames(duplis),1,19)%in%x])))
#identify samples with barcode 
colnames(temp)=j
colnames(prefi)=substr(colnames(prefi),1,19)
#joint matrices
final=cbind(prefi,temp)
write.table(final,"methyM.tsv",sep='\t',quote=F)
#2511M when tar.gz

###########check final data#################################
library(NOISeq)

#u're dragging batch effects, use them as covariates when DM
noiseqData = NOISeq::readData(data = final, factor=subtypeLUAD)
myPCA = NOISeq::dat(noiseqData, type = "PCA", norm = T, logtransf = T)
#logtransf=F fails & norm=F flattens points
png("PCA_methy_global.png")
print({explo.plot(myPCA, samples = c(1,2), plottype = "scores",
                  factor = "subtype")})
dev.off()
png("PCA_methy_global_stage.png")
print({explo.plot(myPCA, samples = c(1,2), plottype = "scores",
                  factor = "tissue")})
dev.off()
### check if we need race and gender? i dont think so

png("PCA_methy_global_race.png")
print({explo.plot(myPCA, samples = c(1,2), plottype = "scores",
                  factor = "race")})
dev.off()

print({explo.plot(myPCA, samples = c(1,2), plottype = "scores",
                  factor = "gender")})