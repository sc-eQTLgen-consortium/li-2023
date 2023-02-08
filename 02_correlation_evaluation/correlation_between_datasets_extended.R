###############################################################################
# In order to better interpret the correlation levels:
# check if correlation between single cell and bulk (ImmuNexUT) is higher for matched
# cell types compared to not matched cell types
###############################################################################

library(data.table)
library(reticulate) # to read the single cell data (numpy)
library(ggplot2)
library(viridis)
library(dplyr)

theme_set(theme_bw())

#Rename cell types
ct_fullname<-setNames(c("CD8+ T cells","monocytes","NK cells","B cells","DC"),
                      c("CD8T","monocyte","NK","B","DC"))

#Path to different single cell dataset
datasets<-c(mio_v3="co-expression_indivs_combined/",
            mio_v2="co-expression_indivs_combined/one_million_version2/",
            stemi_v2="co-expression_indivs_combined/stemi/version2/",
            stemi_v3="co-expression_indivs_combined/stemi/version3/",
            pilot="co-expression_indivs_combined/ng_updated_version/")

#File endings for different single cell datasets
file_suffixes<-c(mio_v3="_UT_correlation.csv",
                 mio_v2="_UT_correlation.csv",
                 stemi_v2="_t8w_correlation.csv",
                 stemi_v3="_t8w_correlation.csv",
                 pilot="_correlation.csv")

#Name on plots for different single cell datasets
dataset_names<-c(mio_v3="Oelen (v3)",
                 mio_v2="Oelen (v2)",
                 stemi_v2="van Blokland (v2)",
                 stemi_v3="van Blokland (v3)",
                 pilot="van der Wijst")

#Different bulk datasets
bulk_datasets<-c("BLUEPRINT","BIOS","ImmuNexUT")

resort<-function(corr){
  #Split into two genes
  corr$gene1<-gsub(";.*","",corr$V1)
  corr$gene2<-gsub(".*;","",corr$V1)
  
  #Order them alphabetically
  corr$V1<-ifelse(corr$gene1 < corr$gene2,corr$V1,
                  paste0(corr$gene2,";",corr$gene1))
  corr$gene1<-NULL
  corr$gene2<-NULL
  
  return(corr)
}

################################################################################
# Compare single cell vs ImmuNexUT - all cell types against all cell types
################################################################################

# Cell type matching
ct_mapping<-data.frame(sc_ct=c("CD4T","CD8T","B","monocyte","NK","DC"),
                       imn_ct=c("Naive_CD4","Naive_CD8","Naive_B","CL_Mono","NK","mDC"))

corr_comp<-NULL
for(i in 1:nrow(ct_mapping)){
  
  ct <- ct_mapping$imn_ct[i]
  #cell_type<- "CD4T"
  
  #Load ImmuNexUT data
  combat_tmm<-fread(paste0("imd_paper_rna_data/norm_count/",ct,"_norm_count.txt"))
  
  #Load the different single cell data sets
  for(c1 in 1:length(datasets)){
    
    #Load for each single cell dataset all cell types
    for(cell_type in ct_mapping$sc_ct){
      
      #Load single cell data set
      corr_c2<-fread(paste0(datasets[c1],cell_type,
                            "/",cell_type,file_suffixes[c1]))
      corr_c2<-resort(corr_c2)
      
      #Filter the ImmuNexUT data set
      expressed_genes<-union(gsub(";.*","",corr_c2$V1),
                             gsub(".*;","",corr_c2$V1))
      expressed_genes<-intersect(expressed_genes,combat_tmm$V1)
      
      combat_tmm_filtered<-combat_tmm[combat_tmm$V1 %in% expressed_genes,]
      combat_tmm_filtered<-as.data.frame(combat_tmm_filtered)
      rownames(combat_tmm_filtered)<-combat_tmm_filtered$V1
      combat_tmm_filtered$V1<-NULL
      
      #Calculation correlation
      cor_matrix<-cor(t(combat_tmm_filtered),method="spearman")
      cor_matrix<-reshape2::melt(cor_matrix)
      cor_matrix$Var1<-as.character(cor_matrix$Var1)
      cor_matrix$Var2<-as.character(cor_matrix$Var2)
      cor_matrix<-cor_matrix[cor_matrix$Var1 < cor_matrix$Var2,]
      cor_matrix$V1<-paste0(cor_matrix$Var1,";",cor_matrix$Var2)
      cor_matrix$Var1<-NULL
      cor_matrix$Var2<-NULL
      
      #Compare BIOS with single cell
      corr<-merge(corr_c2,cor_matrix,by=c("V1"))
      
      #Unique genes
      num_genes<-length(union(gsub(";.*","",corr$V1),
                              gsub(".*;","",corr$V1)))
      
      corr_comp<-rbind(corr_comp,
                       data.frame(sc_ct=cell_type,
                                  bulk_ct=ct,
                                  c1=dataset_names[c1],
                                  c2="ImmuNexUT",
                                  gene_pairs=nrow(corr),
                                  genes_unique=num_genes,
                                  corr=cor(corr[[2]],corr[[3]],method="pearson")))
    }
  }
}

# Save correlations
write.table(corr_comp,
            file="co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_immunexut_mixedcts.tsv",
            sep="\t",row.names = FALSE,quote=FALSE)


################################################################################
# Plot the results
################################################################################

corr_comp<-fread("co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_immunexut_mixedcts.tsv")

ct_mapping<-data.frame(sc_ct=c("CD4T","CD8T","B","monocyte","NK","DC"),
                       imn_ct=c("Naive_CD4","Naive_CD8","Naive_B","CL_Mono","NK","mDC"))

#Order single cell and bulk the same way
corr_comp$sc_ct<-factor(corr_comp$sc_ct,levels=ct_mapping$sc_ct)
corr_comp$bulk_ct<-factor(corr_comp$bulk_ct,levels=ct_mapping$imn_ct)

g<-ggplot(corr_comp,aes(x=sc_ct,y=bulk_ct,fill=corr))+
  geom_tile()+
  geom_text(aes(label=paste0(round(corr,3),"\n(",genes_unique,")"),
                color=ifelse(corr<0.6,'white','black')),size=3)+
  scale_color_manual(values=c("black","white"))+
  facet_wrap(~c1)+
  xlab("Cell type - single cell")+
  ylab("Cell type - bulk")+
  scale_fill_viridis("Correlation",limits=c(0,1))+
  theme(legend.position = "bottom")+
  guides(color="none")
  

print(g)

ggsave(g,file="correlation_mixed_cts.png",height=7,width=9)

################################################################################
# Normalize the columns to always by the diagonal (matched cell types)
################################################################################

ct_mapping_list<-setNames(c("CD4T","CD8T","B","monocyte","NK","DC"),
                          c("Naive_CD4","Naive_CD8","Naive_B","CL_Mono","NK","mDC"))
corr_comp$bulk_matched_ct<-ct_mapping_list[corr_comp$bulk_ct]
corr_diagonal<-corr_comp[corr_comp$sc_ct==corr_comp$bulk_matched_ct,c("sc_ct","c1","corr")]
colnames(corr_diagonal)<-c("sc_ct","c1","diag_corr")

corr_comp<-merge(corr_comp,corr_diagonal,by=c("sc_ct","c1"))
corr_comp$rel_corr<-corr_comp$corr/corr_comp$diag_corr

g<-ggplot(corr_comp,aes(x=sc_ct,y=bulk_ct,fill=rel_corr))+
  geom_tile()+
  geom_text(aes(label=paste0(round(rel_corr,3),"\n(",round(corr,3),")"),
                color=ifelse(rel_corr<1,'white','black')),size=3)+
  scale_color_manual(values=c("black","white"))+
  facet_wrap(~c1)+
  xlab("Cell type - single cell")+
  ylab("Cell type - bulk")+
  scale_fill_viridis("Relative corr")+
  theme(legend.position = "bottom")+
  guides(color="none")


print(g)

ggsave(g,file="correlation_mixed_cts_normalized.png",height=7,width=9)
