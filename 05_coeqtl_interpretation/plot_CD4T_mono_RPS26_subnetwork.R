# ------------------------------------------------------------------------------
# Create network of co-eQTLs connected with rs1131017-RPS26 in CD4+ T cells
# and/or Monocytes,  color edges by direction of effect,
# annotate them with GO terms for transcription initiation
# and lympocyte activity and transfer it to cytoscape for final layouting
# Input: coeQTL results from CD4+ T cells and Monocytes
# Output: cytoscape graph
# ------------------------------------------------------------------------------

library(data.table)
library(igraph) #to plot the graph
library(dplyr)
library(ggplot2)
library(AnnotationHub) # for GO annotation
library(GO.db) # for GO annotation
library(RCy3) # to move igraph object to cytoscape

theme_set(theme_bw())

# Load current set of coeQTL
coeqtls_CD4T<-fread("coeqtl_mapping/CD4T/coeqtls_fullresults_fixed_withAF.sig.tsv")
coeqtls_mono<-fread("coeqtl_mapping/monocyte/coeqtls_fullresults_fixed_withAF.sig.tsv")

#Correct for the direction of effect
coeqtls_CD4T$MetaPZ<-ifelse(coeqtls_CD4T$AF>0.5,(-1)*coeqtls_CD4T$MetaPZ,
                            coeqtls_CD4T$MetaPZ)
coeqtls_mono$MetaPZ<-ifelse(coeqtls_mono$AF>0.5,(-1)*coeqtls_mono$MetaPZ,
                            coeqtls_mono$MetaPZ)

coeqtls<-rbind(coeqtls_CD4T,coeqtls_mono)

#Get gene1 and gene2 correctly sorted
coeqtls$gene1<-gsub(";.*","",coeqtls$Gene)
coeqtls$gene2<-gsub(".*;","",coeqtls$Gene)
coeqtls$eqtlgene<-gsub(".*_","",coeqtls$snp_eqtlgene)
coeqtls$gene2<-ifelse(coeqtls$gene1 == coeqtls$eqtlgen, coeqtls$gene2,
                      coeqtls$gene1)
coeqtls$gene1<-coeqtls$eqtlgene
coeqtls$direction<-ifelse(coeqtls$MetaPZ>0,"positive","negative")

#Set the direction is NA for the ones without matching direction
coeqtls<-unique(coeqtls[,c("snp_genepair","snp_eqtlgene","gene2","direction")])
non_matching_dir<-intersect(coeqtls$snp_genepair[coeqtls$direction=="positive"],
                            coeqtls$snp_genepair[coeqtls$direction=="negative"])
coeqtls$direction[coeqtls$snp_genepair %in% non_matching_dir]<-"not_maching"

#Remove duplicate entries for non-matching directions
coeqtls<-unique(coeqtls[,c("snp_genepair","snp_eqtlgene","gene2","direction")])
coeqtls$type<-ifelse(coeqtls$snp_genepair %in% coeqtls_CD4T$snp_genepair,
                     ifelse(coeqtls$snp_genepair %in% coeqtls_mono$snp_genepair,
                            "both","CD4T"),"mono")
# Filter for RPS26
coeqtls<-coeqtls[coeqtls$snp_eqtlgene=="rs1131017_RPS26",]

#Check how much the correlation structure matches for overlapping examples
corr_RPS26_cd4t<-fread("coeqtl_interpretation/correlation_structure_coeqtl_rps26.tsv")
corr_RPS26_mono<-fread("coeqtl_interpretation/correlation_structure_coeqtl_rps26monocyte.tsv")
corr_comp<-merge(corr_RPS26_cd4t,corr_RPS26_mono,by=c("Var1","Var2"))

g<-ggplot(corr_comp,aes(x=value.x,y=value.y))+
  geom_point()+geom_abline()+
  xlab("Correlation CD4T")+ylab("Correlation Monocytes")
print(g)
ggsave("coeqtl_interpretation/plots_filtered/correlation_structure_RPS26_compare_cts.png")

#Load interaction structure
corr_RPS26<-rbind(corr_RPS26_cd4t,corr_RPS26_mono)
corr_RPS26<-corr_RPS26%>%
  group_by(Var1,Var2)%>%
  summarise(value=max(value))%>%
  as.data.frame()

#First thing, show only strong associations
corr_RPS26<-corr_RPS26[corr_RPS26$value>0.2,]
corr_RPS26$value<-NULL
colnames(corr_RPS26)<-c("snp_eqtlgene","gene2")
corr_RPS26$direction<-"correlation"

#Edges
# edges_combined<-rbind(coeqtls[,c("snp_eqtlgene","gene2","direction")],
#                   corr_RPS26)
edges_combined<-coeqtls[,c("snp_eqtlgene","gene2","direction")]

# Graph object
graph_object <- graph_from_edgelist(as.matrix(edges_combined[,1:2]),
                                    directed=FALSE)

#Color SNPs and genes differently
V(graph_object)$node_type<-ifelse(startsWith(names(V(graph_object)),"rs"),
                                  "eQTL","gene2")

#Get for nodes also to which type of network they belong
nodes_cd4t<-unique(unlist(coeqtls[coeqtls$type %in% c("CD4T","both"),
                                  c("snp_eqtlgene","gene2")]))
nodes_mono<-unique(unlist(coeqtls[coeqtls$type %in% c("mono","both"),
                                  c("snp_eqtlgene","gene2")]))
nodes_both<-intersect(nodes_cd4t,nodes_mono)

V(graph_object)$coeqtl<-ifelse(names(V(graph_object)) %in% nodes_both,
                               "both",
                               ifelse(names(V(graph_object)) %in% nodes_cd4t,
                                      "CD4T","mono"))

#Color edges
E(graph_object)$direction<-edges_combined$direction

# Check how many are associated with transcription initiation (GO:0006413)
# and T cell activiation (GO:0042110) / lympocyte activity (GO:0046649)
# or any of its offspring

ah <- AnnotationHub()
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]

go_transcriptinit <- AnnotationDbi::select(orgdb, c("GO:0006413",GOBPOFFSPRING[["GO:0006413"]]), 
                                  "SYMBOL", "GO")
go_tcell <- AnnotationDbi::select(orgdb, c("GO:0046649",GOBPOFFSPRING[["GO:0046649"]]), 
                                  "SYMBOL", "GO")

#intersect(go_transcriptinit$hgnc_symbol,go_tcell$hgnc_symbol)
V(graph_object)$GO<-ifelse(names(V(graph_object)) %in% go_transcriptinit$SYMBOL,
                           ifelse(names(V(graph_object)) %in% go_tcell$SYMBOL,
                                  "both","transcript_init"),
                           ifelse(names(V(graph_object)) %in% go_tcell$SYMBOL,
                                  "tcell_active","none"))

coeqtls$enrichment<-ifelse(coeqtls$gene2 %in% go_transcriptinit$SYMBOL,
                           ifelse(coeqtls$gene2 %in% go_tcell$SYMBOL,
                                  "both","transcript_init"),
                           ifelse(coeqtls$gene2 %in% go_tcell$SYMBOL,
                                  "tcell_active","none"))
table(coeqtls$enrichment,coeqtls$direction)

plot(graph_object,vertex.label=NA, vertex.size=3)

# Import into cytoscape
cytoscapePing()
createNetworkFromIgraph(graph_object,"RPS26_network_function")

