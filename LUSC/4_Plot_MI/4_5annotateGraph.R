setwd("~/lungsquamouscells/MI/filteredCys/")
library(RCy3)
library(tidyverse)
########################PARAMETERS & PACKAGES
#net=commandArgs(trailingOnly=TRUE)
net="hsa04512.LUSC.cys"
# id=unlist(strsplit(net,'.',fixed=T))[1]
#subty=unlist(strsplit(net,'.',fixed=T))[2]
id="hsa04512"
subty="LUSC"

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(igraph))
library(RCy3)

#########################LOAD NET
openSession(net)
g=createIgraphFromNetwork(subty)
names=data.frame(cbind("name"=V(g)$name,"label"=V(g)$label))
#######################################3HIGHLIGHT FUNCTIONAL NODES
#get the cluster that contains the function
# chosen=read_tsv("../subsamples/chosen.tsv",show_col_types=F)
chosen_bp=read_tsv("/home/mdiaz/lungsquamouscells/functional_enrichment/bp.enrichment",show_col_types=F)
chosen_kegg=read_tsv("/home/mdiaz/lungsquamouscells/functional_enrichment/kegg.enrichment", show_col_types=F)
chosen=rbind(chosen_bp, chosen_kegg)

chosen=chosen%>%filter(ID==id&subtype==subty)
cl=chosen$group
funs=chosen$Description

#get related functions
groups=read_tsv("/home/mdiaz/lungsquamouscells/functional_enrichment/Groups_per_component.tsv",show_col_types=F)
groups=groups%>%filter(subtype==subty,group==cl)
if(nrow(groups)>1){
  funs=groups%>%dplyr::select(Description)%>%unlist}

#get "functional" nodes
print("Identify functional nodes")
setwd("/home/mdiaz/lungsquamouscells/functional_enrichment/")
enriched=list(bp=read_tsv("bp.enrichment",show_col_types=F),
              kegg=read_tsv("kegg.enrichment",show_col_types=F))
setwd("~/lungsquamouscells/MI/filteredCys/")
known_genes=lapply(enriched,function(x) 
  x%>%dplyr::filter(Description%in%funs&subtype==subty)%>%
    dplyr::select(Description,geneID)%>%separate_rows(geneID,sep='/'))
#get symbols for kegg
if(nrow(known_genes$kegg)>0){
  library(biomaRt)
  mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",
                  version=105)
  myannot <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id',
                                'entrezgene_id'),mart=mart)
  colnames(myannot)[3]="geneID"
  known_genes$kegg=merge(known_genes$kegg,myannot,by="geneID",all.x=T)%>%
    dplyr::select(Description,hgnc_symbol)
  colnames(known_genes$kegg)[2]="geneID"
  known_genes=do.call(rbind,known_genes)
}else{
  known_genes=known_genes$bp	
}

known_genes=known_genes%>%filter(geneID%in%unique(unlist(strsplit(V(g)$label,","))))
funs=funs[funs%in%known_genes$Description]
#don't use labels to merge or it'll mess ids
colnames(known_genes)[2]="label"
names= names %>% separate_rows(label, sep=",")
known_genes=merge(known_genes,names,by="label")
gk=graph.data.frame(known_genes[,2:3])#funs & names
#fix labels
i=known_genes%>%distinct(label,name)
i=rbind(i,cbind(label=funs,name=funs))
i2= duplicated(i$name)
i2= i[!i2, ]
i=i2
V(gk)$label=i$label[order(match(i$name,V(gk)$name))]
print("Functions net")
try(createNetworkFromIgraph(gk,"known"))
#add function annotation
print("Annotating functional nodes")
try(mergeNetworks(sources = paste("known", subty, sep = ","),title="merged"))
print("Nice network")

#function - nodes edges are a different kind of edge
i=selectNodes(nodes=funs,by.col="name",network="merged")$nodes
#node and edge names are necesary to set bypass
setNodeColorBypass(node.names=i,new.colors='#EBEBEB',network="merged")
#repeated or it will die
selectNodes(nodes=funs,by.col="name",network="merged")
# $nodes
# [1] 2343
# 
# $edges
# [1] 2351 2349
i=selectEdgesAdjacentToSelectedNodes(network="merged")$edges
setEdgeLineStyleBypass(edge.names=i,new.styles="EQUAL_DASH")

#################################################HIGHLIGHT TFS
print("checking TFs")
tfs <- read_csv("/home/mdiaz/lungsquamouscells/functional_enrichment/DatabaseExtract_v_1.01.csv")
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
  mirIDs=read_tsv("/home/mdiaz/lungsquamouscells/functional_enrichment/miR_IDs.tsv",skip=0)
  mirIDs=mirIDs[mirIDs$precursor%in%query,]
  miRtargets=get_multimir(mirna=c(query,mirIDs$mature),
                          summary=F,table="predicted",legacy.out=F)
  miRtargets=multiMiR::select(miRtargets,keys="predicted",
                              columns=columns(miRtargets),keytype="type")
  miRtargets[miRtargets$target_symbol%in%names$label,]
}
