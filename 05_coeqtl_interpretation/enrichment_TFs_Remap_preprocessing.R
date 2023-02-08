# ------------------------------------------------------------------------------
# Preparation of enrichment of TFBS among co-eGenes using Remap annotations:
# filter Remap file for blood related cell lines 
# (due to size of file this takes some time)
# Input: bed file from Remap 2022 (from https://remap.univ-amu.fr/download_page)
# Output: bed file filtered for blood cell lines 
# ------------------------------------------------------------------------------

library(rtracklayer)

#Read peak file from Remap2022
peaks<-import("tfbs_enrichment_remap/remap2022_nr_macs2_hg19_v1_0.bed.gz")

# Filter blood related cell lines
ann <- t(matrix(unlist(strsplit(values(peaks)[,"name"], ":", fixed=T)), nrow=2))
colnames(ann) <- c("TF", "conditions")
ann <- as.data.frame(ann,stringsAsFactors=FALSE)

conditions<-data.frame(condition=unique(unlist(strsplit(ann$conditions,","))),
                       blood_related=FALSE)
conditions<-conditions[order(conditions$condition),]

blood_related_terms<-c("ALL","AML",
                       "B-cell","BJAB","BL41","blood",
                       "CLL","DC","erythroid",
                       "erythroid-progenitor","GM",
                       "Jurkat","K-562","Kasumi",
                       "LCL","leukemia","lymphoblast","lymphocyte",
                       "macrophage","MM1-S","monocyte",
                       "neutrophil","P493","peripheral-blood",
                       "SEM","T-cell","Th1","Th17","THP-1","U-937",
                       "monocyte")

for(term in blood_related_terms){
  conditions[grep(term, conditions$condition),"blood_related"] <- TRUE
}

write.table(conditions,file="tfbs_enrichment_remap/conditions_remap2022.tsv",
            quote=FALSE,sep="\t")

conditions<-conditions[conditions$blood_related,]

ann$blood_related<-FALSE
for(term in blood_related_terms){
  ann[grep(term, ann$conditions),"blood_related"] <- TRUE
}

values(peaks)<-ann
peaks<-peaks[peaks$blood_related]
peaks$blood_related<-NULL

#Put back into name column as only this column is exported from rtracklayer
peaks$name<-paste0(peaks$TF,":",peaks$conditions)
export(peaks,"tfbs_enrichment_remap/remap2022_nr_macs2_hg19_v1_0_blood_related.bed")
