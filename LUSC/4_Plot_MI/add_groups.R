setwd("/home/mdiaz/lungsquamouscells/functional_enrichment")

library(tidyverse)
groups <- readRDS("groups.RDS")

bp <- read_tsv("BP.enrichment")
kegg <- read_tsv("KEGG.enrichment")
names(groups) <- unique(bp$subtype)

add_cluster <- function(id, subtype, comp, type){
  if(type == "bp"){
    description <- bp$Description[bp$ID == id & bp$subtype == subtype & bp$component == comp]
  } else if(type == "kegg"){
    description <- kegg$Description[kegg$ID == id & kegg$subtype == subtype & kegg$component == comp]
  } else {
    cat("type incorrect")
  }
  
  # select group based on subtype 
  g1 <- as.data.frame(groups[subtype])
  colnames(g1) <- "description"
  return(g1[description,])
}

add_clusters <- Vectorize(add_cluster)
bp$group <- add_clusters(bp$ID, bp$subtype, bp$component, "bp")
kegg$group <- add_clusters(kegg$ID,kegg$subtype, kegg$component, "kegg")

#write_tsv(bp, "BP.enrichment")
#write_tsv(kegg, "KEGG.enrichment")

write_tsv(bp, "bp.enrichment")
write_tsv(kegg, "kegg.enrichment")