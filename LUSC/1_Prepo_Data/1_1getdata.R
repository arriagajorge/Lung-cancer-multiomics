setwd("/home/mdiaz/lungsquamouscells/prepo_data/")

library(BiocManager)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(VennDiagram)
library(data.table)
library(tidyverse)
library(dplyr)

##########SAMPLE IDs PER DATA TYPE#####################
mthyltnLUSC <-  GDCquery(project = "TCGA-LUSC",
                         data.category = "DNA Methylation",
                         platform="Illumina Human Methylation 450")

df1 <- mthyltnLUSC[[1]][[1]]
mthyltnLUSC = getResults(mthyltnLUSC)

i=substr(mthyltnLUSC$cases,1,19)

xprssnLUSC <- GDCquery(project = "TCGA-LUSC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type="STAR - Counts")#
xprssnLUSC=getResults(xprssnLUSC)

j=substr(xprssnLUSC$cases,1,19)
mirnasLUSC <- GDCquery(project = "TCGA-LUSC",
                       data.category = "Transcriptome Profiling",
                       data.type = "miRNA Expression Quantification")
mirnasLUSC=getResults(mirnasLUSC)
k=substr(mirnasLUSC$cases,1,19)

##############CONCOURRENT MEASURES########################
sapply(list(i,j,k),function(x) length(unique(x)))
#[1] 412 553 523 
#only samples with concurrent measurements are useful
samples=intersect(intersect(i,j),k)
length(samples)
#[1] 368    

samples_df=data.frame(cbind(samples,
                            sapply(samples,function(x) 
                              unique(as.character(mthyltnLUSC$sample_type[i==x])))))
colnames(samples_df)[2]="tissue"
# unique(samples_df$tissue)
# [1] "Primary Tumor"       "Solid Tissue Normal"
samples_df$patient=substr(samples_df$samples,1,12) 

samples <- samples_df

subtypes=TCGAquery_subtype(tumor="LUSC")#subtype per patient
sum(samples$patient%in%subtypes$patient)
# [1] 75
# unique(subtypes$Expression.Subtype)
# [1] basal     classical secretory primitive

#only classified samples  are useful
samples=samples[samples$patient%in%subtypes$patient,]
table(subtypes$Expression.Subtype[subtypes$patient%in%samples$patient])
# basal classical primitive secretory 
# 13        27        10        22 

#normal subtype is discarded NOT normal tissue 
samples_2<-samples[which(!samples$patient%in%subtypes$patient[   #samples_2 is only a name for dont overwrite
  subtypes$Expression.Subtype=="Normal"]|
    samples$tissue=="Solid Tissue Normal"),] #there not are Subtype Normal

temp=table(samples_2[samples_2$patient%in%samples_2$patient[duplicated(
  samples_2$patient)],2:3])
table(apply(temp,2,paste,collapse=""))

#subtype per sample
samples_2$subtype=sapply(samples_2$patient,function(x) 
  subtypes$Expression.Subtype[subtypes$patient==x])

levels(samples_2$subtype) <- c(levels(samples_2$subtype), "normal")
samples_2 <- samples_2 %>% 
  mutate(subtype = replace(subtype, tissue =="Solid Tissue Normal", "normal"))

table(samples_2$subtype)
# basal classical primitive secretory    normal 
# 13        27        10        22         3 
#there are few information per subtype
#delete samples from subtype normal
samples_2 <- samples_2 %>% 
  mutate(subtype = ifelse(subtype %in% c("basal", "classical", "secretory", "primitive"), "LUSC", subtype),
         subtype = ifelse(tissue == "Solid Tissue Normal", "normal", subtype))


write.table(samples_2,"subtypeLUSC.tsv",sep='\t',quote=F,row.names=F)

#getwd()
#plot intersections
i=i[mthyltnLUSC$sample_type!="Solid Tissue Normal"]
j=j[xprssnLUSC$sample_type!="Solid Tissue Normal"]
k=k[mirnasLUSC$sample_type!="Solid Tissue Normal"]
levels_lung <- as.character(unique(subtypes$expression_subtype))

lista <- lapply(levels_lung,function(x)
  venn.diagram(x=list(
    A=unique(j[substr(j,1,12)%in%subtypes$patient[subtypes$expression_subtype==x]]),
    B=unique(i[substr(i,1,12)%in%subtypes$patient[subtypes$expression_subtype==x]]),
    C=unique(k[substr(k,1,12)%in%subtypes$patient[subtypes$expression_subtype==x]])),
    col = "transparent", fill = c("cornflowerblue","green","red"),
    alpha = 0.50,cex = 1.5,cat.cex = 1.5,margin = 0.1,
    fontfamily = "sans",cat.fontfamily=rep("sans",3),
    category.names=c("RNAseq","HM450","miRNAseq"),filename=x)) #estos nombres estan bien?

venn.diagram(x=list(
  A=substr(mthyltnLUSC$cases[mthyltnLUSC$sample_type=="Solid Tissue Normal"],1,19),
  B=substr(xprssnLUSC$cases[xprssnLUSC$sample_type=="Solid Tissue Normal"],1,19),
  C=substr(mirnasLUSC$cases[mirnasLUSC$sample_type=="Solid Tissue Normal"],1,19)),
  col = "transparent", fill = c("cornflowerblue","green","red"),
  alpha = 0.50,cex = 1.5,cat.cex = 1.5,margin = 0.1,
  fontfamily = "sans",cat.fontfamily=rep("sans",3),
  category.names=c("RNAseq","HM450","miRNAseq"),filename="Normal")


###############ADD CLINICAL INFO##########################
clin <- GDCquery_clinic("TCGA-LUSC","clinical")
#change tumor_stage by ajcc_pathologic_stage 
clin <- clin[,c("bcr_patient_barcode","gender",
                "ajcc_pathologic_stage","race","vital_status")]
samples_2=cbind(samples_2,t(sapply(samples_2$patient,function(x) 
  clin[clin$bcr_patient_barcode==x,2:4])))
table(clin$gender) 
# female   male 
# 131    373 
table(clin$ajcc_pathologic)
# Stage I   Stage IA   Stage IB   Stage II  Stage IIA  Stage IIB  Stage III Stage IIIA Stage IIIB   Stage IV 
# 3         90        152          3         65         95          3         63         19          7 
table(clin$race)
# asian black or african american              not reported                     white 
# 9                        31                       113                       351 

samples_ <- as.matrix(samples_2)
setwd("/home/mdiaz/lungsquamouscells/prepo_data/")
write.table(samples_,"subtypeLUSC.tsv",sep='\t',quote=F,row.names=F)