################################################################################
# Compare effect sizes (Z-scores) calculated in each individual dataset  
# (before the meta-analysis)
# Input: coeqtls results of the respective cell type
# Output: pairwise plot showing the differences for each combination of cohorts
################################################################################

library(GGally) #to generate pairwise comparison plots
library(viridis)

coeqtl_dir<-"coeqtl_mapping/output/filtered_results/"
plot_dir<-"coeqtl_interpretation/plots_filtered/"  

cell_type<-"CD4T"

# Load current set of coeQTL
coeqtls<-fread(paste0(coeqtl_dir,"UT_",
                      cell_type,"/coeqtls_fullresults.all.tsv.gz"))
coeqtls$gene1<-gsub(";.*","",coeqtls$Gene)
coeqtls$gene2<-gsub(".*;","",coeqtls$Gene)

# Gene 1 and 2 should be ordered alphabetically, but there is an issue regarding
# small and capital letters (so order them again!)
coeqtls$swap<-ifelse(coeqtls$gene1 > coeqtls$gene2,coeqtls$gene1,coeqtls$gene2)
coeqtls$gene1<-ifelse(coeqtls$gene1 > coeqtls$gene2,coeqtls$gene2,coeqtls$gene1)
coeqtls$gene2<-coeqtls$swap
coeqtls$swap<-NULL

# Filter for significant coeQTLs
sign_coeqtls<-coeqtls[coeqtls$gene2_isSig == "TRUE" &
                        coeqtls$snp_qval <= 0.05,]

print(paste(nrow(sign_coeqtls),"significant coeQTLs from",
            nrow(coeqtls),"pairs"))
print(paste("CoeQTLs consisting of:",
            length(unique(sign_coeqtls$GeneSymbol)), "unique gene pairs from",
            length(unique(c(sign_coeqtls$gene1,sign_coeqtls$gene2))),"unique genes",
            "and",length(unique(sign_coeqtls$SNP)),"unique SNPs"))

# Check Z score distribution
z_scores<-strsplit(sign_coeqtls$`DatasetZScores(ng;onemillionv2;onemillionv3;stemiv2)`,
                   split=";")
z_scores<-matrix(as.numeric(unlist(z_scores)),ncol=4,byrow=TRUE)
z_scores<-as.data.frame(z_scores)
colnames(z_scores)<-c("ng","onemillionv2","onemillionv3","stemiv2")
z_scores$coeqtl<-sign_coeqtls$snp_genepair
z_scores$meta_z<-sign_coeqtls$MetaPZ

#Flip the Z-scores so that AF is always representing the minor allele
z_scores$AF<-sign_coeqtls$SNPEffectAlleleFreq
for(colN in c("ng","onemillionv2","onemillionv3","stemiv2","meta_z")){
  z_scores[,colN]<-ifelse(z_scores$AF>=0.5,z_scores[,colN]*(-1),z_scores[,colN])
}

#Rename Z score columns
colnames(z_scores)[1:4]<-c("van der Wijst","Oelen (v2)","Oelen (v3)", "van Blokland (v2)")
z_scores<-z_scores[,c("Oelen (v2)","Oelen (v3)", "van Blokland (v2)",
                      "van der Wijst","coeqtl","meta_z","AF")]
#Plot comparison of Z scores between cohorts
lowerfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping)+
    geom_bin2d()+
    scale_fill_viridis("Density",breaks=c(2,7),labels = c("Low", "High"))+
    geom_hline(yintercept=0)+geom_vline(xintercept=0)
}

g<-ggpairs(z_scores[1:4],
           lower=list(continuous=wrap(lowerfun)),
           legend=c(2,1))+
  theme(legend.position = "bottom")
ggsave(g,file=paste0(plot_dir,cell_type,
                     "_zscore_dist_cohorts.pdf"),
       height=6,width=6)