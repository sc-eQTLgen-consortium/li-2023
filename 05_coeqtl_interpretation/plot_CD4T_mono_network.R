# ------------------------------------------------------------------------------
# Create network of all co-eQTLs from CD4+ T cells
# and/or Monocytes (displaying only the large connected component) 
# and color edges by direction of effect
# Input: coeQTL results from CD4+ T cells and Monocytes
# Output: cytoscape graph
# ------------------------------------------------------------------------------

library(data.table)
library(igraph) #to plot the graph
library(RCy3) # to transfer to cytoscape

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

# Graph object
graph_object <- graph_from_edgelist(as.matrix(coeqtls[,c("snp_eqtlgene","gene2")]),
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

V(graph_object)$coeqtl<-ifelse(names(V(graph_object))%in% nodes_both,
                               "both",
                               ifelse(names(V(graph_object)) %in% nodes_cd4t,
                                      "CD4T","mono"))

#Color edges
E(graph_object)$direction<-coeqtls$direction

#Remove all the small components
comps<-components(graph_object)
subgraph<-induced_subgraph(graph_object,vids=names(which(comps$membership==1)))
plot(subgraph,vertex.label=NA, vertex.size=3)

# Transfer results to cytoscape
cytoscapePing()
createNetworkFromIgraph(subgraph,"eqtl_network")

