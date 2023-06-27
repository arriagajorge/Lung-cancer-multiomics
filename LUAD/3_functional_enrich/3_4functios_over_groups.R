#install.packages('ape')
#install.packages('phytools')

library(ape)
library(phytools)
library(tidyverse)
library(igraph)
library(RCy3)


# if(!"RCy3" %in% installed.packages()){
#   install.packages("BiocManager")
#   BiocManager::install("RCy3")
# }

#library(RCy3)
options("cyRestApi"="http://localhost:1234")
cytoscapePing ()

setwd("/home/jvasquez/Documents/Omics/LUAD/3_FuncEnrichment/")
bp=read_tsv("BP.enrichment")
k=read_tsv("KEGG.enrichment")
enriched=list(BP=bp,KEGG=k)

##############FIND "CROSSLINKED" FUNCTIONS
coenriched=do.call(rbind,enriched)%>%group_by(subtype)%>%
  group_map(~table(.x[,c("Description","component")]))
#matrix with the component intersections
intersection=lapply(coenriched,function(x) crossprod(t(x)))
#matrix with the component unions
union=lapply(coenriched,function(z) sapply(1:nrow(z),function(x) 
  sapply(1:nrow(z),function(y) sum(colSums(z[c(x,y),])>0))))
#Jaccard index for the components
coenriched=lapply(1:1,function(x) intersection[[x]]/union[[x]])

#don't forget 1-x to have identical sets together
trees=lapply(coenriched,function(x) hclust(as.dist(1-x)))

#get the groups that are enriched exactly in the same components
groups=lapply(trees,function(x) cutree(x,h=0))
saveRDS(groups, "groups.RDS")
groups=do.call(rbind,lapply(1:1,function(y) 
  data.frame(cbind("subtype"=unique(enriched$BP$subtype)[y],
                   "Description"=names(groups[[y]]),
                   "group"=groups[[y]]))))
write_tsv(groups,"Groups_per_component.tsv")

#copied from http://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/Phylogenetic-trees.nb.html
names(trees) <- c("LUAD")
temp=ape::as.phylo(trees$LUAD)
ig=ape::as.igraph.phylo(temp,FALSE)
ig=set_edge_attr(ig,"distance",value=temp$edge.length)
createNetworkFromIgraph(ig)
createColumnFilter('junctions', 'id', "^Node\\\\d+$", "REGEX")
junctions<-getSelectedNodes()
setNodeWidthBypass(junctions,1)
setNodeHeightBypass(junctions,1)
setNodeLabelBypass(junctions, "")

