#!/usr/bin/env Rscript

########################PARAMETERS & PACKAGES
#net=commandArgs(trailingOnly=TRUE)
setwd("/home/jvasquez/Documents/Omics/LUAD/4_MI/sort/filt/")
# net=commandArgs(trailingOnly=TRUE)
#net = "GO:0000768.LUAD.cys" #no cool #id=cg19003797 problem 
#net = "GO:0014902.LUAD.cys" #cool #id=CAPN2 
#net = "GO:0045666.LUAD.cys" #4 vertex but not known genes #id=miR-199a-2
#net = "GO:0072028.LUAD.cys" #no cool #id=MYOF
#net = "GO:0072210.LUAD.cys" #no cool #id=MYOF
#net = "GO:0001656.LUAD.cys" #no cool #id=MYOF
#net = "GO:0021534.LUAD.cys" #no cool #id=miR-125b-2
#net = "GO:0045682.LUAD.cys" #coool #id=miR-125b-2
#net = "GO:0072074.LUAD.cys" #no cool #id=MYOF
#net = "GO:0072224.LUAD.cys" #no cool #id=MYOF
#net = "GO:0001657.LUAD.cys" #no cool #id=MYOF
#net = "GO:0021885.LUAD.cys" #4 vertex but not known genes #id=miR-199a-2
#net = "GO:0060142.LUAD.cys" #cool #id=CAPN2 
#net = "GO:0072075.LUAD.cys" #no cool #id=MYOF
#net = "GO:0097696.LUAD.cys" #cool id=miR-125b-1-2 
#net = "GO:0001658.LUAD.cys" #no cool #id=MYOF
#net = "GO:0021924.LUAD.cys" #no cool #id RAI2
#net = "GO:0060143.LUAD.cys" #cool #id=CAPN2 
#net = "GO:0072078.LUAD.cys" #no cool #id=MYOF
#net = "GO:0140253.LUAD.cys" #cool #id=CAPN2 
#net = "GO:0001823.LUAD.cys" #no cool #id=MYOF
#net = "GO:0021930.LUAD.cys" #no cool #id RAI2
#net = "GO:0060675.LUAD.cys" #no cool #id=MYOF
#net = "GO:0072080.LUAD.cys" #no cool #id=MYOF
#net = "GO:1901739.LUAD.cys" #cool #id=CAPN2
#net = "GO:0002320.LUAD.cys" #cool #id=miR-125b-2
#net = "GO:0021936.LUAD.cys" #no cool #id=RAI2
#net = "GO:0060993.LUAD.cys" #no cool #id=MYOF
#net = "GO:0072088.LUAD.cys" #no cool #id=MYOF
#net = "GO:1901741.LUAD.cys" #cool id=CAPN2
#net = "GO:0002328.LUAD.cys" #cool #id=miR-125b-2
#net = "GO:0022029.LUAD.cys" #no #miR-199a-2 #4 vertex
#net = "GO:0061005.LUAD.cys" #no cool id=MYOF
#net = "GO:0072163.LUAD.cys" #no cool id=MYOF
#net = "hsa04928.LUAD.cys" #no cool #id=miR-500a
#net = "GO:0006949.LUAD.cys" #cool #id=CAPN2
#net = "GO:0032835.LUAD.cys" #no cool #id=MYOF
#net = "GO:0061326.LUAD.cys" #no cool #id=MYOF
#net = "GO:0072164.LUAD.cys" #no cool #id=MYOF
#net = "hsa05206.LUAD.cys" #cool #id MIR-199A-2
#net = "GO:0007520.LUAD.cys" #cool #id=CAPN2
#net = "GO:0042531.LUAD.cys" #cooool miR-125b-2
#net = "GO:0061333.LUAD.cys" #no cool #id=MYOF
#net = "GO:0072171.LUAD.cys" #no cool #id=MYOF
#net = "GO:0008344.LUAD.cys" #4 vertex but not known genes #id=miR-199a-2
#net = "GO:0045604.LUAD.cys" #coooool #id=MIR125B1
#net = "GO:0072009.LUAD.cys" #no cool #id=MYOF
#net = "GO:0072202.LUAD.cys" #no cool #id=MYOF

id = unlist(strsplit(net,'.',fixed=T))[1]
subty = unlist(strsplit(net,'.',fixed=T))[2]

library(RCy3)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(igraph))

#########################LOAD NET
openSession(net)
g=createIgraphFromNetwork(subty)
names=data.frame(cbind("name"=V(g)$name,"label"=V(g)$label))

#######################################3HIGHLIGHT FUNCTIONAL NODES
#get the cluster that contains the function
# chosen=read_tsv("../subsamples/chosen.tsv",show_col_types=F)
chosen_bp=read_tsv("/home/jvasquez/Documents/Omics/LUAD/3_FuncEnrichment/BP.enrichment",show_col_types=F)
chosen_kegg=read_tsv("/home/jvasquez/Documents/Omics/LUAD/3_FuncEnrichment/KEGG.enrichment",show_col_types=F)
chosen = rbind(chosen_bp, chosen_kegg)

chosen=chosen%>%filter(ID==id&subtype==subty)
cl=unique(chosen$group)
funs=chosen$Description

#get related functions
groups=read_tsv("/home/jvasquez/Documents/Omics/LUAD/3_FuncEnrichment/Groups_per_component.tsv",show_col_types=F)
groups=groups%>%filter(subtype%in%subty,group%in%cl)
if(nrow(groups)>1){
  funs=groups%>%dplyr::select(Description)%>%unlist}

#get "functional" nodes
print("Identify functional nodes")
enriched=list(BP=read_tsv("/home/jvasquez/Documents/Omics/LUAD/3_FuncEnrichment/BP.enrichment",show_col_types=F),
              KEGG=read_tsv("/home/jvasquez/Documents/Omics/LUAD/3_FuncEnrichment/KEGG.enrichment",show_col_types=F))
known_genes=lapply(enriched,function(x) 
  x%>%dplyr::filter(Description%in%funs&subtype==subty)%>%
    dplyr::select(Description,geneID)%>%separate_rows(geneID,sep='/'))
#get symbols for KEGG
if(nrow(known_genes$KEGG)>0){
  library(biomaRt)
  mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
                  version=105)
  myannot <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id',
                                'entrezgene_id'),mart=mart)
  colnames(myannot)[3]="geneID"
  known_genes$KEGG=merge(known_genes$KEGG,myannot,by="geneID",all.x=T)%>%
    dplyr::select(Description,hgnc_symbol)
  colnames(known_genes$KEGG)[2]="geneID"
  known_genes=do.call(rbind,known_genes)
}else{
  known_genes=known_genes$BP	
}

known_genes
V(g)$label

label_fix = toupper(gsub("-", "",  unlist(strsplit(V(g)$label,",")) ) )
V(g)$label = label_fix
known_genes$geneID = toupper(gsub("-", ""  , known_genes$geneID ) )
known_genes=known_genes%>%filter(geneID%in%unique(label_fix))
names$label=toupper(gsub("-","", names$label))

known_genes=known_genes%>%filter(geneID%in%unique(unlist(strsplit(V(g)$label,","))))
if(length(known_genes$geneID) < 1){
  stop("no known genes")
}
funs=funs[funs%in%known_genes$Description]
#don't use labels to merge or it'll mess ids
colnames(known_genes)[2]="label"
known_genes=merge(known_genes,names,by="label")
gk=graph.data.frame(known_genes[,2:3])#funs & names
#fix labels
i=known_genes%>%distinct(label,name)
i=rbind(i,cbind(label=funs,name=funs))
V(gk)$label=unique( i$label[order(match(i$name,V(gk)$name))] )
print("Functions net")
try(createNetworkFromIgraph(gk,"known"))
#add function annotation
print("Annotating functional nodes")
# try(mergeNetworks(c("known", subty),"merged"))
# try(mergeNetworks(sources=c("known", subty),"merged"))
# try(mergeNetworks(sources = list("known", subty), name = "merged"))
try(mergeNetworks(sources = paste("known", subty, sep = ","),title="merged"))
print("Nice network")

#function - nodes edges are a different kind of edge
i=selectNodes(nodes=funs,by.col="name",network="merged")$nodes
#node and edge names are necesary to set bypass
setNodeColorBypass(node.names=i,new.colors='#EBEBEB',network="merged")
#repeated or it will die
selectNodes(nodes=funs,by.col="name",network="merged")
i=selectEdgesAdjacentToSelectedNodes(network="merged")$edges
setEdgeLineStyleBypass(edge.names=i,new.styles="EQUAL_DASH")

#################################################HIGHLIGHT TFS
print("checking TFs")
tfs <- read_csv("/home/jvasquez/Documents/Omics/LUAD/4_MI/sort/filt/DatabaseExtract_v_1.01.csv") #download from http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv
tfs <- tfs[tfs$`Is TF?`=="Yes",]
tfs <- tfs$`HGNC symbol`
#tfs=readLines("humanTFsLambert2018.list")
tfs=names%>%filter(label%in%tfs)%>%distinct(name)%>%unlist
if(length(tfs) > 0){
  setNodeBorderColorBypass(tfs, new.colors='#ab52eb',network="merged")
  setNodeBorderWidthBypass(tfs, 8,network="merged")
}
# setNodeBorderColorBypass(tfs, new.colors='#ab52eb',network="merged")
# setNodeBorderWidthBypass(tfs, 8,network="merged")

saveSession(filename=net)
####################CHECK DBs
print("Getting wikigene description")
#get gene description from biomaRt for all nodes
library(biomaRt)
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
                version=105)
#https://dec2021.archive.ensembl.org
myannot=getBM(attributes = c("hgnc_symbol","wikigene_description"),
              filter="hgnc_symbol",values=names$label,mart=mart)

#PUBMEAD
library(rentrez)

print("Checking functions in pubmed")
#m=subgraph.edges(g,E(g)[from(names$name[names$label%in%known_genes$geneID])|to(names$name[names$label%in%known_genes$geneID])])
#query=names%>%filter(name%in%V(m)$name&!label%in%known_genes$geneID)%>%dplyr::select(label)%>%unlist
#search only the nodes that aren't responsible for the enrichment
query=names$label[!names$label%in%known_genes$label]
query=unique(unlist(strsplit(query,",")))
knownfun=lapply(funs,function(x) #for every function in the net
  lapply(query,function(y) #and every non-functional feature
    #search papers connecting them
    entrez_search(db="pubmed",term=paste(x,y,sep=' AND '))))
write_tsv(as.data.frame(do.call(rbind,lapply(knownfun,function(x) do.call(rbind,lapply(x,function(y) cbind(paste(y$ids,collapse=','),y$QueryTranslation)))))),
          file=paste(subty,"knownfun",sep='.'))

print("Checking edges on pubmed")
edges=data.frame(get.edgelist(g))
colnames(edges)[2]="name"
edges=merge(edges,names,by="name")[,2:3]
colnames(edges)[2]="pair1"
colnames(edges)[1]="name"
edges=unique(merge(edges,names,by="name")[,2:3])
colnames(edges)[2]="pair2"
edges=edges%>%separate_rows(pair2,sep=',',convert=T)%>%separate_rows(pair2,sep=',',convert=T)
known=apply(edges,1,function(x) 
  entrez_search(db="pubmed",
                term=paste(x[1],x[2],sep=' AND ')))
known=known[sapply(known,function(x) x$count)>0]
write_tsv(as.data.frame(do.call(rbind,lapply(known,function(x) cbind(paste(x$ids,collapse=' OR '),x$QueryTranslation)))),
          file=paste(subty,"knownedges",sep='.'))

#multimiR miRNA-target prediction
query=query[grep("miR",query)]
if(length(query)>0){
  library(multiMiR)
  query=gsub("miR","hsa-mir",query)
  mirIDs=read_tsv("miR_IDs.tsv",skip=0)
  mirIDs=mirIDs[mirIDs$precursor%in%query,]
  miRtargets=get_multimir(mirna=c(query,mirIDs$mature),
                          summary=F,table="predicted",legacy.out=F)
  miRtargets=multiMiR::select(miRtargets,keys="predicted",
                              columns=columns(miRtargets),keytype="type")
  miRtargets[miRtargets$target_symbol%in%names$label,]
}

########################WHICH TYPE OF EDGES ARE MORE SHARED
#########################EXCLUSIVE EDGES

#down 1D91C0
#up E31A1C