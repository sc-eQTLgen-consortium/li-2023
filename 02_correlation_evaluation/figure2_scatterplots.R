# ------------------------------------------------------------------------------
# Create inset plots for Main Figure 2 (a,b,d), showing scatterplots of
# gene pair-wise Spearman correlation values between two data sets for
# a) Oelen v3 dataset vs van Blokland v2 dataset (both CD4+ T cells)
# b) ImmuNexUT - van Blokland v2 (naive CD4+ T cells and CD4+ T cells)
# c) Blueprint - ImmuNexUT (both naive CD4+ T cells)
# ------------------------------------------------------------------------------

library(data.table)
library(reticulate) # to read the single cell data (numpy)
library(reshape2)
library(ggplot2)
library(viridis)
library(ggpubr)

np <- import("numpy")

theme_set(theme_bw())

#Load single cell
load_sc_corr_data<-function(path){
  
  # Load single cell data
  corr_ct<-fread(path)
  corr_ct$gene1<-sapply(corr_ct$V1,function(s) strsplit(s,";")[[1]][1])
  corr_ct$gene2<-sapply(corr_ct$V1,function(s) strsplit(s,";")[[1]][2])
  corr_ct$V1<-NULL
  
  #Order so that gene1 is always the one first in alphabet
  corr_ct$swap<-ifelse(corr_ct$gene1 > corr_ct$gene2,corr_ct$gene1,corr_ct$gene2)
  corr_ct$gene1<-ifelse(corr_ct$gene1 > corr_ct$gene2,corr_ct$gene2,corr_ct$gene1)
  corr_ct$gene2<-corr_ct$swap
  corr_ct$swap<-NULL
  
  return(corr_ct)
}

#Load data saved in numpy format
load_numpy_data<-function(path,rowname_suffix,colname_suffix,corr_sc){
  corr_bios <- np$load(paste0(path,"npy"), mmap_mode="r")
  row_names<-fread(paste0(path,rowname_suffix),header=FALSE)
  rownames(corr_bios)<-row_names$V1
  col_names<-fread(paste0(path,colname_suffix),header=FALSE)
  colnames(corr_bios)<-col_names$V1
  rm(row_names,col_names)
  
  #Filter for single cell data
  sc_genes<-sort(union(corr_sc$gene1,corr_sc$gene2))
  sc_genes<-sc_genes[sc_genes %in% colnames(corr_bios)]
  corr_bios<-corr_bios[sc_genes,sc_genes]
  corr_bios<-reshape2::melt(corr_bios)
  corr_bios$Var1<-as.character(corr_bios$Var1)
  corr_bios$Var2<-as.character(corr_bios$Var2)
  colnames(corr_bios)[1:2]<-c("gene1","gene2")
  
  
  #Order so that gene1 is always the one first in alphabet
  corr_bios$swap<-ifelse(corr_bios$gene1 > corr_bios$gene2,corr_bios$gene1,corr_bios$gene2)
  corr_bios$gene1<-ifelse(corr_bios$gene1 > corr_bios$gene2,corr_bios$gene2,corr_bios$gene1)
  corr_bios$gene2<-corr_bios$swap
  corr_bios$swap<-NULL
  corr_bios<-corr_bios[corr_bios$gene1!=corr_bios$gene2,]
  
  return(corr_bios)
}

#Create ggplot 2d histogram based on two correlation data frames
create_corr_plot<-function(corr_d1,corr_d2,
                           xlab_text,ylab_text,
                           annot_text_size=9,annot_text_digits=3){
  
  #Merge both
  corrs<-merge(corr_d1,corr_d2,by=c("gene1","gene2"))
  
  print(paste("Overlapping genes:",length(union(corrs$gene1,
                                                corrs$gene2))))
  
  corr_corr<-cor(corrs$corr1,corrs$corr2,
                 method="pearson")
  
  #Plot
  g<-ggplot(corrs,aes(corr1,corr2))+
    geom_bin2d(bins=50)+
    xlab(xlab_text)+
    ylab(ylab_text)+
    xlim(-1,1)+ylim(-1,1)+
    scale_fill_distiller("Density",palette="BuPu",trans="log10",
                         breaks = c(2, 600), 
                         labels = c("Low", "High"))+
    annotate(geom="text", x=-0.95, y=0.95,size=annot_text_size,
             hjust = 0,vjust=1,
             label=paste0("r = ",format(corr_corr,digits=annot_text_digits)))+
    ggtitle("Pairwise gene correlation")+
    geom_smooth(method="lm",color="black")+
    theme(legend.position="none",
          plot.title=element_text(size=25),
          axis.title=element_text(size=25),
          axis.text=element_text(size=20))
  
  return(g)
}
  
  
################################################################################
# For 2a: Oelen v3 - van Blokland v2
################################################################################

main_celltype<-"CD4T"

corr_oelen<-load_sc_corr_data(paste0("co-expression_indivs_combined/",
                                     main_celltype,"/",main_celltype,
                                     "_UT_correlation.csv"))

corr_stemi<-load_sc_corr_data(paste0("co-expression_indivs_combined/stemi/version2/",
                                     main_celltype,"/",main_celltype,
                                     "_t8w_correlation.csv"))

colnames(corr_oelen)<-c("corr1","gene1","gene2")
colnames(corr_stemi)<-c("corr2","gene1","gene2")

#Create gggplot
g<-create_corr_plot(corr_oelen,corr_stemi,
                    xlab_text="Oelen (v3)", ylab_text="van Blokland (v2)")

g_leg<-get_legend(g+theme(legend.position = "bottom"))
g_leg<-as_ggplot(g_leg)
ggsave(g_leg,file="bios/plots/figure2_legend_inset.pdf",width=3,height=1)


ggsave(g,file="bios/plots/figure2a_exampleplot.pdf",width=5,height=5)

################################################################################
# For 2b: ImmuNexUT - van Blokland v2
################################################################################

corr_immu<-fread("imd_paper_rna_data/correlation/Naive_CD4_correlation.txt")
corr_immu$V1<-NULL

colnames(corr_immu)<-c("gene1","gene2","corr1")

#Create gggplot
g<-create_corr_plot(corr_immu,corr_stemi,
                    xlab_text="ImmuNexUT",ylab_text="van Blokland (v2)")

ggsave(g,file="bios/plots/figure2b_exampleplot.pdf",width=5,height=5)

################################################################################
# For 2c: Blueprint - ImmuNexUT
################################################################################

#Load Blueprint data
corr_bp<-load_numpy_data(path="blueprint_data/tcel_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.spearmanR.",
                         rowname_suffix="rows.txt",
                         colname_suffix="cols.txt",
                         corr_immu)

colnames(corr_bp)<-c("gene1","gene2","corr2")

#Create gggplot
g<-create_corr_plot(corr_immu,corr_bp,
                    xlab_text="ImmuNexUT",ylab_text="BLUEPRINT")
ggsave(g,file="bios/plots/figure2c_exampleplot.pdf",width=5,height=5)
