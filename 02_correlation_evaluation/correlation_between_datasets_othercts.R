# ------------------------------------------------------------------------------
# Extension of correlation_between_datasets.R (which looks only at CD4+ T cells)
# for other cell types: get Pearson correlation between data sets
# * for single cell vs single cell data set
# * for single cell vs bulk data set
# * for bulk vs bulk data set (here only monocytes)
# Plot one heatmap for each comparison
# ------------------------------------------------------------------------------

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
# 1) Compare single cell with each other
################################################################################

corr_comp<-NULL
for(cell_type in c("CD8T","monocyte","NK","B","DC")){
  
  for(c1 in 1:(length(datasets)-1)){
    
    #Read correlation file one
    dataset_name1<-dataset_names[c1]
    corr_c1<-fread(paste0(datasets[c1],cell_type,
                          "/",cell_type,file_suffixes[c1]))
    corr_c1<-resort(corr_c1)
    
    #Unique genes
    num_genes<-length(union(gsub(";.*","",corr_c1$V1),
                            gsub(".*;","",corr_c1$V1)))
    
    corr_comp<-rbind(corr_comp,
                     data.frame(cell_type,
                                c1=dataset_name1,
                                c2=dataset_name1,
                                gene_pairs=nrow(corr_c1),
                                genes_unique=num_genes,
                                corr=1))
    
    for(c2 in (c1+1):length(datasets)){
      
      #Read correlation file two
      dataset_name2<-dataset_names[c2]
      corr_c2<-fread(paste0(datasets[c2],cell_type,"/",
                            cell_type,file_suffixes[c2]))
      corr_c2<-resort(corr_c2)
      
      corr<-merge(corr_c1,corr_c2,by=c("V1"))
      
      #Unique genes
      num_genes<-length(union(gsub(";.*","",corr$V1),
                              gsub(".*;","",corr$V1)))
      
      corr_comp<-rbind(corr_comp,
                       data.frame(cell_type,
                                  c1=dataset_name1,
                                  c2=dataset_name2,
                                  gene_pairs=nrow(corr),
                                  genes_unique=num_genes,
                                  corr=cor(corr[[2]],corr[[3]],method="pearson")))
    }
  }
  
  
  #Read correlation file one
  c1<-length(datasets)
  dataset_name1<-dataset_names[c1]
  corr_c1<-fread(paste0(datasets[c1],cell_type,
                        "/",cell_type,file_suffixes[c1]))
  
  #Unique genes
  num_genes<-length(union(gsub(";.*","",corr_c1$V1),
                          gsub(".*;","",corr_c1$V1)))
  
  corr_comp<-rbind(corr_comp,
                   data.frame(cell_type,
                              c1=dataset_name1,
                              c2=dataset_name1,
                              genes_unique=num_genes,
                              gene_pairs=nrow(corr_c1),
                              corr=1))
}

# Save correlations
write.table(corr_comp,
            file="co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_datasets_othercts.tsv",
            sep="\t",row.names = FALSE,quote=FALSE)

################################################################################
# Plot comparison of single cell vs single cell (Supplementary Figure)
################################################################################

corr_comp<-fread("co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_datasets_othercts.tsv")
corr_comp$c1<-factor(corr_comp$c1,levels=dataset_names)
corr_comp$c2<-factor(corr_comp$c2,levels=dataset_names)
corr_comp$cell_type<-ct_fullname[corr_comp$cell_type]

g<-ggplot(corr_comp,aes(x=c1,y=c2,fill=corr))+
  geom_tile()+
  geom_text(aes(label=paste0(round(corr,3),"\n(",genes_unique,")"),
                color=ifelse(corr<0.6,'white','black')),size=3)+
  scale_color_manual(values=c("black","white"))+
  xlab("Single cell data set")+
  ylab("Single cell data set")+
  scale_fill_viridis("Correlation",limits=c(0,1))+
  facet_wrap(~cell_type)+
  scale_y_discrete(labels=c("Oelen (v3)","Oelen (v2)","van Blokland\n(v2)",
                            "van Blokland\n(v3)","van der Wijst"))+
  scale_x_discrete(labels=c("Oelen\n(v3)","Oelen\n(v2)","van\nBlokland\n(v2)",
                            "van\nBlokland\n(v3)","van der\nWijst"))+
  theme(legend.position=c(0.9,0.1))+
  guides(color=FALSE)
print(g)

ggsave(g,file=paste0("co-expression_indivs_combined/plots/corr_single_cell_othercts.png"),
       width=8.5,height=6.5)

#Get also CD4 T cell results
corr_comp_ct<-fread("co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_datasets.tsv")
corr_comp_ct$cell_type<-"CD4T"
corr_comp_ct<-corr_comp_ct[,colnames(corr_comp),with=FALSE]
corr_comp<-rbind(corr_comp_ct,corr_comp)

#Get median correlation
corr_comp%>%
  group_by(cell_type)%>%
  summarise(mean(corr),median(corr),min(corr),max(corr))

#Overall distribution across all cell types
summary(corr_comp$corr)

################################################################################
# 2) Compare single cell with bulk
################################################################################

################################################################################
# For Blueprint - Monocytes
################################################################################

#Blueprint Monocyte correlation
path<-"blueprint_data/mono_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.1PCAsOverSamplesRemoved.spearmanr."
rowname_suffix<-"genes.txt"

corr_c1 <- np$load(paste0(path,"npy"), mmap_mode="r")
row_names<-fread(paste0(path,rowname_suffix),header=FALSE)
rownames(corr_c1)<-row_names$V1
colnames(corr_c1)<-row_names$V1
rm(row_names)

corr_comp<-NULL
cell_type<-"monocyte"
for(c1 in 1:length(datasets)){
  
  #Load single cell data set
  corr_c2<-fread(paste0(datasets[c1],cell_type,
                        "/",cell_type,file_suffixes[c1]))
  corr_c2<-resort(corr_c2)
  
  #Filter the Blueprint data set
  expressed_genes<-union(gsub(";.*","",corr_c2$V1),
                         gsub(".*;","",corr_c2$V1))
  expressed_genes<-intersect(expressed_genes,colnames(corr_c1))
  
  corr_c1_filtered<-corr_c1[expressed_genes,expressed_genes]
  corr_c1_filtered<-reshape2::melt(corr_c1_filtered)
  corr_c1_filtered$Var1<-as.character(corr_c1_filtered$Var1)
  corr_c1_filtered$Var2<-as.character(corr_c1_filtered$Var2)
  corr_c1_filtered<-corr_c1_filtered[corr_c1_filtered$Var1 < corr_c1_filtered$Var2,]
  corr_c1_filtered$V1<-paste0(corr_c1_filtered$Var1,";",corr_c1_filtered$Var2)
  corr_c1_filtered$Var1<-NULL
  corr_c1_filtered$Var2<-NULL
  
  corr<-merge(corr_c1_filtered,corr_c2,by=c("V1"))
  
  #Unique genes
  num_genes<-length(union(gsub(";.*","",corr$V1),
                          gsub(".*;","",corr$V1)))
  
  corr_comp<-rbind(corr_comp,
                   data.frame(cell_type,
                              c1=dataset_names[c1],
                              c2="BLUEPRINT",
                              gene_pairs=nrow(corr),
                              genes_unique=num_genes,
                              corr=cor(corr[[2]],corr[[3]],method="pearson")))
}

# Save correlations
write.table(corr_comp,
            file="co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_blueprint_mono.tsv",
            sep="\t",row.names = FALSE,quote=FALSE)

################################################################################
# For Blueprint - CD4T
################################################################################

#Blueprint CD4T correlation
path<-"blueprint_data/tcel_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.spearmanR."
rowname_suffix<-"rows.txt"

corr_c1 <- np$load(paste0(path,"npy"), mmap_mode="r")
row_names<-fread(paste0(path,rowname_suffix),header=FALSE)
rownames(corr_c1)<-row_names$V1
colnames(corr_c1)<-row_names$V1
rm(row_names)

corr_comp<-NULL
cell_type<-"CD4T"
for(c1 in 1:length(datasets)){
  
  #Load single cell data set
  corr_c2<-fread(paste0(datasets[c1],cell_type,
                        "/",cell_type,file_suffixes[c1]))
  corr_c2<-resort(corr_c2)
  
  #Filter the Blueprint data set
  expressed_genes<-union(gsub(";.*","",corr_c2$V1),
                         gsub(".*;","",corr_c2$V1))
  expressed_genes<-intersect(expressed_genes,colnames(corr_c1))
  
  corr_c1_filtered<-corr_c1[expressed_genes,expressed_genes]
  corr_c1_filtered<-reshape2::melt(corr_c1_filtered)
  corr_c1_filtered$Var1<-as.character(corr_c1_filtered$Var1)
  corr_c1_filtered$Var2<-as.character(corr_c1_filtered$Var2)
  corr_c1_filtered<-corr_c1_filtered[corr_c1_filtered$Var1 < corr_c1_filtered$Var2,]
  corr_c1_filtered$V1<-paste0(corr_c1_filtered$Var1,";",corr_c1_filtered$Var2)
  corr_c1_filtered$Var1<-NULL
  corr_c1_filtered$Var2<-NULL
  
  corr<-merge(corr_c1_filtered,corr_c2,by=c("V1"))
  
  #Unique genes
  num_genes<-length(union(gsub(";.*","",corr$V1),
                          gsub(".*;","",corr$V1)))
  
  corr_comp<-rbind(corr_comp,
                   data.frame(cell_type,
                              c1=dataset_names[c1],
                              c2="BLUEPRINT",
                              gene_pairs=nrow(corr),
                              genes_unique=num_genes,
                              corr=cor(corr[[2]],corr[[3]],method="pearson")))
}

# Save correlations
write.table(corr_comp,
            file="co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_blueprint_cd4t.tsv",
            sep="\t",row.names = FALSE,quote=FALSE)

################################################################################
# For ImmuNexUT - all cell types
################################################################################

# Cell type matching
ct_mapping<-data.frame(sc_ct=c("CD4T","CD8T","B","monocyte","NK","DC"),
                       imn_ct=c("Naive_CD4","Naive_CD8","Naive_B","CL_Mono","NK","mDC"))

corr_comp<-NULL
for(i in 1:nrow(ct_mapping)){
  
  ct <- ct_mapping$imn_ct[i]
  cell_type<- ct_mapping$sc_ct[i]
  
  #Load ImmuNexUT data
  combat_tmm<-fread(paste0("imd_paper_rna_data/norm_count/",ct,"_norm_count.txt"))
  
  #Load the different single cell data sets
  for(c1 in 1:length(datasets)){
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
                     data.frame(cell_type,
                                c1=dataset_names[c1],
                                c2="ImmuNexUT",
                                gene_pairs=nrow(corr),
                                genes_unique=num_genes,
                                corr=cor(corr[[2]],corr[[3]],method="pearson")))
  }
}

# Save correlations
write.table(corr_comp,
            file="co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_immunexut_allcts.tsv",
            sep="\t",row.names = FALSE,quote=FALSE)

################################################################################
# For BIOS - all cell types
################################################################################

#Load the bios expression matrix
bios_data<-fread("bios/gene_read_counts_BIOS_and_LLD_passQC.tsv.SampleSelection.ProbesWithZeroVarianceRemoved.TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.noLLDeep.scGeneOnly.txt.gz")

corr_comp<-NULL
for(cell_type in c("CD4T","CD8T","monocyte","NK","B","DC")){
  for(c1 in 1:length(datasets)){
    
    corr_c1<-fread(paste0(datasets[c1],cell_type,
                          "/",cell_type,file_suffixes[c1]))
    corr_c1<-resort(corr_c1)
    expressed_genes<-union(gsub(";.*","",corr_c1$V1),
                            gsub(".*;","",corr_c1$V1))
    
    #Filter BIOS for the genes expressed in the respective data set  
    bios_data_filtered<-bios_data[bios_data$genename %in% expressed_genes,]
    bios_data_filtered<-as.data.frame(bios_data_filtered)
    rownames(bios_data_filtered)<-bios_data_filtered$genename
    bios_data_filtered$genename<-NULL
    
    #Calculate correlation for BIOS
    cor_matrix<-cor(t(bios_data_filtered),method="spearman")
    cor_matrix<-reshape2::melt(cor_matrix)
    cor_matrix$Var1<-as.character(cor_matrix$Var1)
    cor_matrix$Var2<-as.character(cor_matrix$Var2)
    cor_matrix<-cor_matrix[cor_matrix$Var1 < cor_matrix$Var2,]
    cor_matrix$V1<-paste0(cor_matrix$Var1,";",cor_matrix$Var2)
    cor_matrix$Var1<-NULL
    cor_matrix$Var2<-NULL
    
    #Compare BIOS with single cell
    corr<-merge(corr_c1,cor_matrix,by=c("V1"))
    
    #Unique genes
    num_genes<-length(union(gsub(";.*","",corr$V1),
                            gsub(".*;","",corr$V1)))
    
    corr_comp<-rbind(corr_comp,
                     data.frame(cell_type,
                                c1=dataset_names[c1],
                                c2="BIOS",
                                gene_pairs=nrow(corr),
                                genes_unique=num_genes,
                                corr=cor(corr[[2]],corr[[3]],method="pearson")))

  }
}

# Save correlations
write.table(corr_comp,
            file="co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_bios_allcts.tsv",
            sep="\t",row.names = FALSE,quote=FALSE)

################################################################################
# Plot comparison of single cell vs bulk (Supplementary Figure)
################################################################################

#Load the different data sets
corr_comp<-rbind(fread("co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_bios_allcts.tsv"),
                 fread("co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_immunexut_allcts.tsv"),
                 fread("co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_blueprint_mono.tsv"),
                 fread("co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_blueprint_cd4t.tsv"))

corr_comp%>%
  group_by(cell_type,c2)%>%
  summarise(mean(corr),median(corr),min(corr),max(corr))
 
#Remove CD4T cells (already shown in the main figure)
corr_comp<-corr_comp[corr_comp$cell_type != "CD4T",]

corr_comp$c1<-factor(corr_comp$c1,levels=dataset_names)

corr_comp$cell_type<-ct_fullname[corr_comp$cell_type]

g<-ggplot(corr_comp,aes(x=c2,y=c1,fill=corr))+
  geom_tile()+
  geom_text(aes(label=paste0(round(corr,3),"\n(",genes_unique,")"),
                color=ifelse(corr<0.6,'white','black')),size=3)+
  scale_color_manual(values=c("black","white"))+
  xlab("Bulk cell data set")+
  ylab("Single cell data set")+
  scale_fill_viridis("Correlation",limits=c(0,1))+
  facet_wrap(~cell_type)+
  scale_y_discrete(labels=c("Oelen (v3)","Oelen (v2)","van Blokland\n(v2)",
                            "van Blokland\n(v3)","van der Wijst"))+
  theme(legend.position=c(0.9,0.1))+
  guides(color=FALSE)
print(g)

ggsave(g,file=paste0("co-expression_indivs_combined/plots/corr_singlevsbulk_othercts.png"),
       width=8.5,height=6.5)

################################################################################
# 3) Compare bulk vs bulk for Monocytes
################################################################################

#Special function to read bulk data as they are not all saved in the same file type
read_bulk_data<-function(dataset_name){
  
  if(dataset_name=="BLUEPRINT"){
    #Blueprint Monocyte correlation
    path<-"blueprint_data/mono_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.1PCAsOverSamplesRemoved.spearmanr."
    rowname_suffix<-"genes.txt"
    
    corr_c1 <- np$load(paste0(path,"npy"), mmap_mode="r")
    row_names<-fread(paste0(path,rowname_suffix),header=FALSE)
    rownames(corr_c1)<-row_names$V1
    colnames(corr_c1)<-row_names$V1
    rm(row_names)
    
    #Filter for single cell data
    ct_single_cell<-"monocyte"
    corr_sc<-fread(paste0("co-expression_indivs_combined/",ct_single_cell,"/",
                          ct_single_cell,"_UT_correlation.csv"))
    corr_sc$gene1<-gsub(";.*","",corr_sc$V1)
    corr_sc$gene2<-gsub(".*;","",corr_sc$V1)
    sc_genes<-union(corr_sc$gene1,corr_sc$gene2)
    sc_genes<-sc_genes[sc_genes %in% colnames(corr_c1)]
    
    corr_c1<-corr_c1[sc_genes,sc_genes]
    corr_c1<-reshape2::melt(corr_c1)
    corr_c1$Var1<-as.character(corr_c1$Var1)
    corr_c1$Var2<-as.character(corr_c1$Var2)
    colnames(corr_c1)[1:3]<-c("gene1","gene2","corr")
    
    #Order so that gene1 is always the one first in alphabet
    corr_c1<-corr_c1[corr_c1$gene1!=corr_c1$gene2,]
    corr_c1$V1<-paste0(corr_c1$gene1,";",corr_c1$gene2)
    corr_c1$gene1<-NULL
    corr_c1$gene2<-NULL
    
    corr_c1<-corr_c1[,c("V1","corr")]
    
  } else if (dataset_name=="BIOS"){

    #Load the bios expression matrix
    bios_data<-fread("bios/gene_read_counts_BIOS_and_LLD_passQC.tsv.SampleSelection.ProbesWithZeroVarianceRemoved.TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.noLLDeep.scGeneOnly.txt.gz")
    
    # Read single cell data to filter for the expressed genes
    cell_type <- "monocyte"
    corr_sc<-fread(paste0("co-expression_indivs_combined/",cell_type,
                          "/",cell_type,"_UT_correlation.csv"))
    corr_sc<-resort(corr_sc)
    expressed_genes<-union(gsub(";.*","",corr_sc$V1),
                           gsub(".*;","",corr_sc$V1))
    
    #Filter BIOS for the genes expressed in the respective data set  
    bios_data_filtered<-bios_data[bios_data$genename %in% expressed_genes,]
    bios_data_filtered<-as.data.frame(bios_data_filtered)
    rownames(bios_data_filtered)<-bios_data_filtered$genename
    bios_data_filtered$genename<-NULL
    
    #Calculate correlation for BIOS
    corr_c1<-cor(t(bios_data_filtered),method="spearman")
    corr_c1<-reshape2::melt(corr_c1)
    corr_c1$Var1<-as.character(corr_c1$Var1)
    corr_c1$Var2<-as.character(corr_c1$Var2)
    corr_c1<-corr_c1[corr_c1$Var1 < corr_c1$Var2,]
    corr_c1$V1<-paste0(corr_c1$Var1,";",corr_c1$Var2)
    corr_c1$Var1<-NULL
    corr_c1$Var2<-NULL
    
  } else { #ImmuNexUT
    corr_c1<-fread("imd_paper_rna_data/correlation/CL_Mono_correlation.txt")
    #all(corr_c1$gene1 < corr_c1$gene2)
    corr_c1$V1<-paste0(corr_c1$gene1,";",corr_c1$gene2)
    corr_c1$gene1<-NULL
    corr_c1$gene2<-NULL
  }
  
  return(corr_c1)
}

#Compare each bulk dataset against all other
corr_comp<-NULL
for(c1 in 1:(length(bulk_datasets)-1)){
  
  #Read correlation file one
  dataset_name1<-bulk_datasets[c1]
  corr_c1<-read_bulk_data(dataset_name1)
  
  #Unique genes
  num_genes<-length(union(gsub(";.*","",corr_c1$V1),
                          gsub(".*;","",corr_c1$V1)))
  
  corr_comp<-rbind(corr_comp,
                   data.frame(c1=dataset_name1,
                              c2=dataset_name1,
                              gene_pairs=nrow(corr_c1),
                              genes_unique=num_genes,
                              corr=1))
  
  for(c2 in (c1+1):length(bulk_datasets)){
    
    #Read correlation file two
    dataset_name2<-bulk_datasets[c2]
    corr_c2<-read_bulk_data(dataset_name2)
    
    corr<-merge(corr_c1,corr_c2,by=c("V1"))
    
    #Unique genes
    num_genes<-length(union(gsub(";.*","",corr$V1),
                            gsub(".*;","",corr$V1)))
    
    corr_comp<-rbind(corr_comp,
                     data.frame(c1=dataset_name1,
                                c2=dataset_name2,
                                gene_pairs=nrow(corr),
                                genes_unique=num_genes,
                                corr=cor(corr[[2]],corr[[3]],method="pearson")))
  }
}

#Read correlation file one
c1<-length(bulk_datasets)
dataset_name1<-bulk_datasets[c1]
corr_c1<-read_bulk_data(dataset_name1)

#Unique genes
num_genes<-length(union(gsub(";.*","",corr_c1$V1),
                        gsub(".*;","",corr_c1$V1)))

corr_comp<-rbind(corr_comp,
                 data.frame(c1=dataset_name1,
                            c2=dataset_name1,
                            gene_pairs=nrow(corr_c1),
                            genes_unique=num_genes,
                            corr=1))

# Save correlations
write.table(corr_comp,
            file="co-expression_indivs_combined/dataset_comp_summary/correlation_bulk_datasets_monocytes.tsv",
            sep="\t",row.names = FALSE,quote=FALSE)

# Save plot
corr_comp$c1<-factor(corr_comp$c1,levels=bulk_datasets)
corr_comp$c2<-factor(corr_comp$c2,levels=bulk_datasets)

g<-ggplot(corr_comp,aes(x=c1,y=c2,fill=corr))+
  geom_tile()+
  geom_text(aes(label=paste0(round(corr,3),"\n(",genes_unique,")"),
                color=ifelse(corr<0.6,'white','black')),size=3)+
  scale_color_manual(values=c("black","white"))+
  xlab("Bulk data set")+
  ylab("Bulk data set")+
  scale_fill_viridis("Correlation",limits=c(0,1))+
  guides(color=FALSE)

ggsave(g,file="co-expression_indivs_combined/plots/corr_bulk_mono.png",
       width=5,height=3)