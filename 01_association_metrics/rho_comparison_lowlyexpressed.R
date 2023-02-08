# ------------------------------------------------------------------------------
# Explore the differences between Spearman correlation and Rho propensity
# for very lowly expressed genes; due to computational very demanding 
# calculation of Rho values for many genes, the comparison is done on a subset
# of 50 very lowly expressed genes (expressed in >0% and <5% of cells) 
# and 50 very highly expressed genes (expressed in >95% of the cells)
# Input: Seurat object of Oelen v3 dataset (UT monocytes)
# Output: scatterplot for comparison
# ------------------------------------------------------------------------------

library(Seurat)
library(propr)
library(Matrix)
library(ggplot2)

theme_set(theme_bw())

#Load complete seurat object
seurat<-readRDS("seurat_objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds")

#Filter for monocytes and UT timepoint
seurat<-seurat[,seurat$cell_type_lowerres == "monocyte"]
seurat<-seurat[,seurat$timepoint == "UT"]

# Get non-zero-ratio of each gene
nozeroratio<-rowMeans(seurat$RNA@data>0)

full_count_matrix<-t(as.matrix(seurat$RNA@data))

#Select 50 very lowly expressed genes and 50 very highly expressed genes
set.seed(1)
low_genes<-sample(names(nozeroratio)[nozeroratio > 0 & nozeroratio < 0.05],50)
high_genes<-sample(names(nozeroratio)[nozeroratio > 0.9],50)

#Calculate rho values
res<-propr::perb(full_count_matrix,
                 select=c(low_genes,high_genes))@matrix

propr<-reshape2::melt(res)
propr$Var1<-as.character(propr$Var1)
propr$Var2<-as.character(propr$Var2)
propr<-propr[propr$Var1 < propr$Var2,]

#Compare with spearman values
spearman<-cor(full_count_matrix[,c(low_genes,high_genes)],method="spearman")
spearman<-reshape2::melt(spearman)
spearman$Var1<-as.character(spearman$Var1)
spearman$Var2<-as.character(spearman$Var2)
spearman<-spearman[spearman$Var1 < spearman$Var2,]

#Combine both into one plot
all(propr$Var1 == spearman$Var1)
all(propr$Var2 == spearman$Var2)

propr$corr<-spearman$value

propr$type<-ifelse(propr$Var1 %in% low_genes,
                   ifelse(propr$Var2 %in% high_genes,"mixed","both_low"),
                   ifelse(propr$Var2 %in% high_genes,"both_high","mixed"))

g<-ggplot(propr,aes(x=corr,y=value,color=type))+
  geom_point(alpha=0.5)+
  xlab("Spearman correlation")+
  ylab("Rho proportionality")+
  xlim(-0.2,1)+ylim(-0.2,1)+
  scale_color_discrete("Expression gene pair")+
  geom_abline()
ggsave(g,file="test_rho.pdf")
