setwd("/home/mdiaz/lungsquamouscells/functional_enrichment")
library(tidyverse)
library(UpSetR)

bp=read_tsv("BP.enrichment",show_col_types = FALSE)
k=read_tsv("KEGG.enrichment",show_col_types = FALSE)

enriched=list(BP=bp,
              KEGG=k)
get_sets=function(enriched_table,exclusive){
  
  #matrix subtypes vs function
  g=table(unique(enriched_table[,c("ID","subtype")]))
  #to get exclusive functions for another time
  if(exclusive==T){
    rowSums(g)
    g=g[rowSums(g==0)==3,]
  }
  #upset() needs a list of IDs
  sets=apply(g,2,function(y) names(which(y>0)))
  return(sets)}

functions=lapply(enriched,get_sets,exclusive=F)

###################EXCLUSIVE FUNCTIONS###########3#################
KEGGenrich <- k
exclusive=lapply(enriched,get_sets,exclusive=F)
write_tsv(unique(KEGGenrich[KEGGenrich$ID%in%unlist(exclusive$KEGG),c("ID","Description")]),"KEGG.exclusive")
# #class comes from https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext&filedir=
# kegg=readLines("br08901.keg")
# kegg=kegg[grep("^[A-Z]",kegg,perl=T)]
# kegg=kegg[substr(kegg,1,1)!="A"]
# i=which(substr(kegg,1,1)=="B")
# i=c(i,length(kegg)+1)
# keggL=lapply(1:59,function(x) kegg[(i[x]+1):(i[x+1]-1)])
# names(keggL)=substr(kegg[i[1:59]],4,nchar(kegg[i[1:59]]))
# keggL=lapply(keggL,function(x) substr(x,13,nchar(x)))
# keggDF=as.data.frame(do.call(rbind,lapply(1:length(keggL),function(x) cbind(names(keggL)[x],keggL[[x]]))))
# colnames(keggDF)=c("class","Description")
# ids=read_tsv("KEGG.exclusive")
# merge(ids,keggDF, by="Description",all.x=T)
# mergeData <- merge(ids,keggDF, by="Description",all.x=T) 
#also simplified classes with ':' to the first term

#data frame it
ids=read_tsv("KEGG.exclusive")
#KEGG.classes=as.data.frame(do.call(rbind,lapply(1:1,function(x) 
#  cbind(names(exclusive$KEGG)[x],exclusive$KEGG[[x]]))))
#colnames(KEGG.classes)=c("subtype","ID")
#KEGG.classes=merge(KEGG.classes,ids,by="ID")
ids$LUSC <- "LUSC"
KEGG.classes <- ids

#get BP categories
library(GSEABase)
library(GO.db)
# as in https://support.bioconductor.org/p/128407/
#and https://support.bioconductor.org/p/83375/
fl="http://current.geneontology.org/ontology/subsets/goslim_agr.obo"
#subset used for humans in PMC6800510
slim <- getOBOCollection(fl)#53 ids only
df = select(GO.db, keys(GO.db), "ONTOLOGY")#ontology of all GOids
table(df$ONTOLOGY[df$GOID%in%ids(slim)])#found all slim ids
#BP CC MF 
#21 16 16 
gomap=as.list(GOBPOFFSPRING)#descendents off every goid
found=names(gomap)[names(gomap)%in%ids(slim)]
# [1] "GO:0000003" "GO:0002376" "GO:0005975" "GO:0006259" "GO:0006629" "GO:0007049" "GO:0007610" "GO:0008283" "GO:0009056"
# [10] "GO:0012501" "GO:0016043" "GO:0016070" "GO:0019538" "GO:0023052" "GO:0030154" "GO:0032502" "GO:0042592" "GO:0050877"
# [19] "GO:0050896" "GO:0051234" "GO:1901135"
sum(found%in%df$GOID[df$ONTOLOGY=="BP"])
#[1] 21 #actually only descendents of BP goids
gomap=gomap[names(gomap)%in%ids(slim)]
#format to easy data frames
slim=as.data.frame(do.call(rbind,lapply(1:21,function(x) 
  cbind(names(gomap)[x],gomap[[x]]))))
colnames(slim)=c("parent","child")
slimnames=as.data.frame(sapply(unique(slim$parent),function(x) 
  Term(GOTERM[[x]])))
slimnames$parent=rownames(slimnames)
colnames(slimnames)[1]="name"
names(exclusive$BP) <- exclusive$BP
names(exclusive$KEGG) <- exclusive$KEGG
#count categories per subtype
BP.classes=as.data.frame(do.call(rbind,lapply(1:1,
                                              function(x) cbind(names(exclusive$BP)[x],exclusive$BP[[x]]))))
colnames(BP.classes)=c("subtype","child")
BP.classes=merge(merge(BP.classes,slim,by="child"),
                 slimnames,by="parent")

############TEST OVER-REPRESENTATION OF EXCLUSIVE FUNCTIONS
library(ggplot2)
bias=function(classes,subtype,class){
  totals=table(classes[,c(class,subtype)])
  s=colSums(totals)
  ps=p.adjust(apply(totals,1,function(x) 
    fisher.test(rbind(x,s-x),simulate.p.value=T)$p.val),"bonferroni")
  #TRUE needed coz > 2 categories
  return(rownames(totals)[ps<0.05])}
bias(KEGG.classes,1,2)
#character(0)
i=bias(BP.classes,1,2)
#character(0)

#ni vale la pena dibujarlo
#png("KEGGexclusive.png")
#KEGG.classes%>%count(subtype,class)%>%
# ggplot(aes(x=n,y=class,fill=subtype))+
# geom_bar(stat="identity",position="fill")+
# annotate("text",x=1.05,y=sort(unique(KEGG.classes$class)),
# 	label=KEGG.classes%>%count(class)%>%dplyr::select(n)%>%unlist)+
# scale_x_continuous(labels=scales::percent)+
# theme(text=element_text(size=18),axis.ticks=element_blank(),
# 	panel.background=element_blank(),legend.title=element_blank())+
# xlab("")+ylab("")+scale_fill_viridis_d(option = "plasma")
#dev.off()

png("BPexclusive.png") 
BP.classes%>%count(subtype,name)%>%
  ggplot(aes(x=n,y=name,fill=subtype))+
  geom_bar(stat="identity",position="fill")+
  scale_x_continuous(labels=scales::percent)+
  theme(text=element_text(size=18),axis.ticks=element_blank(),
        panel.background=element_blank(),legend.title=element_blank(),
        legend.position="bottom",legend.margin=margin(-25,0,0,-150))+
  xlab("")+ylab("")+scale_fill_viridis_d(option = "plasma")+
  annotate("text",x=1.08,y=sort(unique(BP.classes$name)),
           label=BP.classes%>%count(name)%>%dplyr::select(n)%>%unlist)+
  annotate("text",y=i,x=-.05,label="*",size=8,vjust=.8)
dev.off()
