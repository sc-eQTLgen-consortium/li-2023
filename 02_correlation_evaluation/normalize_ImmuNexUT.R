# ------------------------------------------------------------------------------
# Normalize ImmuNexUT data (separately for each cell type with a matching
# single-cell cell type) following the description in the corresponding
# publication (filtering lowly expressed genes, TMM normalization and 
# batch correction)
# followed by correlation calculation for all genes expressed in 50% of the cells
# of the Oelen v3 dataset (for comparison with single cell data)
# Input: Count matrices downloaded from 
#        https://humandbs.biosciencedbc.jp/en/hum0214-v5#E-GEAD-397,
#        correlation estimates from Oelen v3 to identify the expressed genes
#        for downstream comparisons
# Output: normalized count matrices (one per cell type), orrelation matrices
#         for all genes expressed in 50% of the cells of the Oelen v3 dataset and
#         plots for comparison between ImmuNexUT and Oelen v3 dataset
# ------------------------------------------------------------------------------

library(data.table)
library(edgeR) #for normalization
library(sva) # for batch correction with combat
library(corrplot) # to plot sample correlations
library(ggplot2)
library(viridis)

theme_set(theme_bw())

# Cell type matching
ct_mapping<-data.frame(sc_ct=c("CD4T","CD8T","B","monocyte","NK","DC"),
                       imn_ct=c("Naive_CD4","Naive_CD8","Naive_B","CL_Mono","NK","mDC"))


for(i in 1:nrow(ct_mapping)){
  
  ct <- ct_mapping$imn_ct[i]
  ct_single_cell<- ct_mapping$sc_ct[i]
  
  # Try to prevent redoing the whole normalization when the correlation
  # is already calculated
  corr_file_name<-paste0("imd_paper_rna_data/correlation/",ct,"_correlation.txt")
  if(! file.exists(corr_file_name)){
    
    counts<-fread(paste0("imd_paper_rna_data/count/",ct,"_count.txt"))
    
    #Format to matrix
    gene_id<-counts$Gene_id
    gene_name<-counts$Gene_name
    counts$Gene_name<-NULL
    counts$Gene_id<-NULL
    counts<-as.matrix(counts)
    rownames(counts)<-gene_name
      
    #Filter lowly expressed genes (at least 10 in > 90% of samples)
    counts<-counts[!(rowSums(counts<10) > 0.9 * ncol(counts)),]
    
    #Normalize using edgeR (TMM plus log-transformed CPM)
    dge <- DGEList(counts=counts)
    dge <- calcNormFactors(dge, method = "TMM")
    tmm <- cpm(dge) #in publication it says log-transformed CPM, but log-transformation
                    #is not working in combination with combat ...
    
    #Remove batch data using combat
    batch_data<-fread("imd_paper_rna_data/clinical_diagnosis_age_sex_v2.txt")
    
    #Filter batch data for samples in the matrix
    batch_data<-batch_data[batch_data$id %in% colnames(tmm),]
    print(paste("Sorted the batch data correctly:",all(batch_data$id == colnames(tmm))))
    
    modcombat = model.matrix(~1, data=batch_data)
    combat_tmm = ComBat(dat=tmm, batch=batch_data$Phase,mod=modcombat, prior.plots = FALSE)
    
    #Check that correlation between samples is high
    cor_matrix<-cor(combat_tmm)
    
    #Filter samples with a correlation coefficient less than 0.9
    cor_coef_mean<-rowMeans(cor_matrix)
    combat_tmm<-combat_tmm[,names(cor_coef_mean)[cor_coef_mean>=0.9]]
    
    # #Plot remaining samples
    # cor_matrix<-cor(combat_tmm)
    # png("imd_paper_rna_data/plots/sample_correlation.png")
    # corrplot(cor_matrix,method="color",order="hclust",tl.col="black",tl.cex=0.2)
    # dev.off()
    
    #Combine genes that appear multiple times in the matrix 
    combat_tmm<- apply(combat_tmm, 2, tapply, rownames(combat_tmm),
                       mean, na.rm=T)
    
    #Save normalized matrix
    write.table(combat_tmm, file=paste0("imd_paper_rna_data/norm_count/",ct,"_norm_count.txt"),
                quote=FALSE,sep="\t")
    
    #Read single cell correlation
    corr_ct<-fread(paste0("co-expression_indivs_combined/",ct_single_cell,"/",
                          ct_single_cell,"_UT_correlation.csv"))
    corr_ct$gene1<-gsub(";.*","",corr_ct$V1)
    corr_ct$gene2<-gsub(".*;","",corr_ct$V1)
    
    #Order so that gene1 is always the one first in alphabet
    corr_ct$swap<-ifelse(corr_ct$gene1 > corr_ct$gene2,corr_ct$gene1,corr_ct$gene2)
    corr_ct$gene1<-ifelse(corr_ct$gene1 > corr_ct$gene2,corr_ct$gene2,corr_ct$gene1)
    corr_ct$gene2<-corr_ct$swap
    corr_ct$swap<-NULL
    corr_ct$V1<-NULL
    
    #Filter for correlation values in CD4 T cells
    expressed_genes<-union(corr_ct$gene1,corr_ct$gene2)
    combat_tmm<-combat_tmm[rownames(combat_tmm) %in% expressed_genes,]
    
    #Calculation correlation
    cor_matrix<-cor(t(combat_tmm),method="spearman")
    cor_matrix<-reshape2::melt(cor_matrix)
    cor_matrix$Var1<-as.character(cor_matrix$Var1)
    cor_matrix$Var2<-as.character(cor_matrix$Var2)
    cor_matrix<-cor_matrix[cor_matrix$Var1 < cor_matrix$Var2,]
    colnames(cor_matrix)<-c("gene1","gene2","corr")
    
    #Save correlation
    write.table(cor_matrix, file=corr_file_name,
                quote=FALSE,sep="\t")
  }
    
  #Compare with single cell correlation
  cor_matrix<-merge(cor_matrix,corr_ct,by=c("gene1","gene2"))
  
  ylab_text<-paste("Correlation ImmuNexUT -",ct)
  plot_path<-paste0("imd_paper_rna_data/plots/correlation_",ct,".png")
  corr_corr<-cor(cor_matrix$UT,cor_matrix$corr)
  
  g<-ggplot(cor_matrix,aes(UT,corr))+
    geom_bin2d(bins=50)+
    xlab("Correlation single cell")+
    ylab(ylab_text)+
    xlim(-1,1)+ylim(-1,1)+
    scale_fill_viridis(trans="log10")+
    annotate(geom="text", x=-0.95, y=0.95,size=8,
             hjust = 0,vjust=1,
             label=paste0("r = ",format(corr_corr,digits=2)))+
    theme(axis.title=element_text(size=16),
          axis.text=element_text(size=14),
          legend.position="none")
  
  ggsave(g,file=plot_path)
}