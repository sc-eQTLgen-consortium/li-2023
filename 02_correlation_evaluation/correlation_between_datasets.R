# ------------------------------------------------------------------------------
# Check Pearson correlation between data sets (for CD4+ T cells)
# * for single cell vs single cell data set
# * for single cell vs bulk data set
# * for bulk vs bulk data set
# Combine all results in one large heatmap
#  -----------------------------------------------------------------------------

library(data.table)
library(reticulate) # to read the single cell data (numpy)
library(ggplot2)
library(viridis)
library(ggpubr)

np <- import("numpy")

theme_set(theme_bw())

cell_type<-"CD4T"

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

bulk_datasets<-c("Blueprint","BIOS","ImmuNexUT")

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
# Compare single cell with each other
################################################################################

corr_comp<-NULL
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
                   data.frame(c1=dataset_name1,
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
                     data.frame(c1=dataset_name1,
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
                 data.frame(c1=dataset_name1,
                            c2=dataset_name1,
                            genes_unique=num_genes,
                            gene_pairs=nrow(corr_c1),
                            corr=1))

corr_comp$c1<-factor(corr_comp$c1,levels=dataset_names)
corr_comp$c2<-factor(corr_comp$c2,levels=dataset_names)

# Save correlations
write.table(corr_comp,
            file="co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_datasets.tsv",
            sep="\t",row.names = FALSE,quote=FALSE)

################################################################################
# Compare single cell with bulk
################################################################################


#Special function to read bulk data as they are not all saved in the same file type
read_bulk_data<-function(dataset_name){
  
  if(dataset_name=="Blueprint"){
    path<-"blueprint_data/tcel_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.spearmanR."
    rowname_suffix<-"rows.txt"
    colname_suffix<-"cols.txt"
      
    corr_c1 <- np$load(paste0(path,"npy"), mmap_mode="r")
    row_names<-fread(paste0(path,rowname_suffix),header=FALSE)
    rownames(corr_c1)<-row_names$V1
    col_names<-fread(paste0(path,colname_suffix),header=FALSE)
    colnames(corr_c1)<-col_names$V1
    rm(row_names,col_names)
    
    #Filter for single cell data
    ct_single_cell<-"CD4T"
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
    corr_c1<-fread("bios/bios_correlation_tcellfiltered.tsv")
    #all(corr_c1$gene1 < corr_c1$gene2)
    corr_c1$V1<-paste0(corr_c1$gene1,";",corr_c1$gene2)
    corr_c1$gene1<-NULL
    corr_c1$gene2<-NULL
  } else { #ImmuNexUT
    corr_c1<-fread("imd_paper_rna_data/correlation/Naive_CD4_correlation.txt")
    #all(corr_c1$gene1 < corr_c1$gene2)
    corr_c1$V1<-paste0(corr_c1$gene1,";",corr_c1$gene2)
    corr_c1$gene1<-NULL
    corr_c1$gene2<-NULL
  }
  
  return(corr_c1)
}

corr_comp<-NULL
for(c1 in 1:length(bulk_datasets)){
  
  dataset_name1<-bulk_datasets[c1]
  corr_c1<-read_bulk_data(dataset_name1)
  
  for(c2 in 1:length(datasets)){
    
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
                     data.frame(c1=dataset_name1,
                                c2=dataset_name2,
                                gene_pairs=nrow(corr),
                                genes_unique=num_genes,
                                corr=cor(corr[[2]],corr[[3]],method="pearson")))
  }
}

# Save correlations
write.table(corr_comp,
            file="co-expression_indivs_combined/dataset_comp_summary/correlation_singlevsbulk_datasets.tsv",
            sep="\t",row.names = FALSE,quote=FALSE)

################################################################################
# Compare bulk with bulk
################################################################################

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
            file="co-expression_indivs_combined/dataset_comp_summary/correlation_bulk_datasets.tsv",
            sep="\t",row.names = FALSE,quote=FALSE)


################################################################################
# Combine all results in one large plot
################################################################################
  
corr_comp<-fread("co-expression_indivs_combined/dataset_comp_summary/correlation_singlecell_datasets.tsv")
corr_comp$c1<-factor(corr_comp$c1,levels=dataset_names)
corr_comp$c2<-factor(corr_comp$c2,levels=dataset_names)
g.1<-ggplot(corr_comp,aes(x=c1,y=c2,fill=corr))+
  geom_tile()+
  geom_text(aes(label=paste0(round(corr,3),"\n(",genes_unique,")")),size=3)+
  xlab("Single cell data set")+
  ylab("Single cell data set")+
  scale_fill_viridis("Correlation",limits=c(0,1))+
  scale_y_discrete(labels=c("Oelen (v3)","Oelen (v2)","van Blokland\n(v2)",
                            "van Blokland\n(v3)","van der Wijst"))+
  scale_x_discrete(labels=c("Oelen\n(v3)","Oelen\n(v2)","van\nBlokland\n(v2)",
                            "van\nBlokland\n(v3)","van der\nWijst"))

corr_comp<-fread("co-expression_indivs_combined/dataset_comp_summary/correlation_singlevsbulk_datasets.tsv")
corr_comp$c1[corr_comp$c1=="Blueprint"]<-"BLUEPRINT"
corr_comp$c2<-factor(corr_comp$c2,levels=bulk_datasets)
corr_comp$c1<-factor(corr_comp$c1,levels=dataset_names)
g.2<-ggplot(corr_comp,aes(x=c2,y=c1,fill=corr))+
  geom_tile()+
  geom_text(aes(label=paste0(round(corr,3),"\n(",genes_unique,")")),size=3,
                color="white")+
  xlab("Bulk data set")+
  ylab("Single cell data set")+
  scale_fill_viridis("Correlation",limits=c(0,1))+
  scale_y_discrete(labels=c("Oelen (v3)","Oelen (v2)","van Blokland\n(v2)",
                            "van Blokland\n(v3)","van der Wijst"))

corr_comp<-fread("co-expression_indivs_combined/dataset_comp_summary/correlation_bulk_datasets.tsv")
corr_comp$c1[corr_comp$c1=="Blueprint"]<-"BLUEPRINT"
corr_comp$c2[corr_comp$c2=="Blueprint"]<-"BLUEPRINT"
corr_comp$c1<-factor(corr_comp$c1,levels=bulk_datasets)
corr_comp$c2<-factor(corr_comp$c2,levels=bulk_datasets)
g.3<-ggplot(corr_comp,aes(x=c1,y=c2,fill=corr))+
  geom_tile()+
  geom_text(aes(label=paste0(round(corr,3),"\n(",genes_unique,")"),
                color=ifelse(corr<0.6,'white','black')),size=3)+
  scale_color_manual(values=c("black","white"))+
  xlab("Bulk data set")+
  ylab("Bulk data set")+
  scale_fill_viridis("Correlation",limits=c(0,1))+
  coord_flip()
  

g_empty<-ggplot()+theme_void()
g<-ggarrange(g.1,g.2,g_empty,g.3,ncol=2,nrow=2,widths=c(4,3),heights=c(4,3),
             common.legend = TRUE,legend="bottom",align="hv")
ggsave(g,file=paste0("co-expression_indivs_combined/plots/corr_datasets_combined.pdf"),
       width=6.5,height=6.5)
