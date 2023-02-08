# ------------------------------------------------------------------------------
# Supplementary figure to show MetaCell overview
# * Expression distribution of genes in a cell
# * number of (meta)cells per sample
# * comparison with Blueprint
# ------------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

theme_set(theme_bw())

# ------------------------------------------------------------------------------
# Expression distribution of genes in a cell
# ------------------------------------------------------------------------------

#Load the single cell object and get expressed genes 
seurat<-readRDS("seurat_objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds")

#Filter for monocytes
seurat<-seurat[,seurat$cell_type_lowerres=="monocyte" &
                 seurat$timepoint=="UT"]

exprGenes.singleCell<-rowSums(as.matrix(seurat@assays$SCT@counts)>0)/ncol(seurat)

print(paste("Number of genes expressed in at least 50% of cells:",
            sum(exprGenes.singleCell>=0.5)))

exprGene.df<-data.frame(expr.perc=sort(exprGenes.singleCell),
                        position=1:length(exprGenes.singleCell),
                        Type="SingleCell",
                        stringsAsFactors = FALSE)

#Metacell
metacell.allsamples<-readRDS("metacell_general/metacell/metacells_K20_minCells10/pseudobulk_metacell.RDS")
#Keep only UT metacells
meta_annot<-read.table("metacell_general/metacell/metacells_K20_minCells10/annotations_metacell.tsv")
#all(meta_annot$metacell == colnames(metacell.allsamples))
metacell.allsamples<-metacell.allsamples[,meta_annot$timepoint == "UT"]
exprGenes.metacell<-rowSums(metacell.allsamples>0)/ncol(metacell.allsamples)

print(paste("Number of genes expressed in at least 50% of cells:",
            sum(exprGenes.metacell>=0.5)))

exprGene.df<-rbind(exprGene.df,
                   data.frame(expr.perc=sort(exprGenes.metacell),
                              position=1:length(exprGenes.metacell),
                              Type="MetaCell",
                              stringsAsFactors = FALSE))

#Leiden
metacell.allsamples<-readRDS("metacell_general/leiden_metacells/metacell_leiden_SCT.RDS")
#Keep only UT metacells
meta_annot<-read.table("metacell_general/leiden_metacells/annotations_mc_leiden_SCT_tp.tsv")
#all(meta_annot$metacell == colnames(metacell.allsamples))
metacell.allsamples<-metacell.allsamples[,meta_annot$condition == "UT"]
exprGenes.metacell<-rowSums(metacell.allsamples>0)/ncol(metacell.allsamples)

print(paste("Number of genes expressed in at least 50% of cells:",
            sum(exprGenes.metacell>=0.5)))

exprGene.df<-rbind(exprGene.df,
                   data.frame(expr.perc=sort(exprGenes.metacell),
                              position=1:length(exprGenes.metacell),
                              Type="Leiden",
                              stringsAsFactors = FALSE))


g.1<-ggplot(exprGene.df,aes(x=position,y=expr.perc,color=Type))+geom_point()+
  xlab("Gene index")+ylab("Expressed in x% of the cells")+
  scale_color_discrete("Method")+
  theme(axis.title = element_text(size=14),
        axis.text=element_text(size=13),
        legend.position = "none")

# ------------------------------------------------------------------------------
# Number of (meta)cells per sample
# ------------------------------------------------------------------------------

counts_all_mc<-NULL

#Load metacell annotation leiden
mc_method<-read.table("metacell_general/leiden_metacells/annotations_mc_leiden_SCT_tp.tsv",
                      stringsAsFactors = FALSE)

mc_method<-mc_method%>%
  group_by(sample,condition)%>%
  summarize(counts=n())

mc_method$method<-"Leiden"
counts_all_mc<-rbind(counts_all_mc,mc_method)

#Load metacell annotation 
mc_method<-read.table("metacell_general/metacell/metacells_K20_minCells10/annotations_metacell.tsv",
                      stringsAsFactors=FALSE)
mc_method<-mc_method%>%
  group_by(sample,timepoint)%>%
  summarize(counts=n())

mc_method$method<-"MetaCell"
colnames(mc_method)<-colnames(counts_all_mc)
counts_all_mc<-rbind(counts_all_mc,mc_method)

#Get number of single cells per sample and condition
sc_annot<-seurat@meta.data
sc_annot<-sc_annot%>%
  group_by(assignment,timepoint)%>%
  summarize(counts=n())
sc_annot$method<-"SingleCell"
colnames(sc_annot)<-colnames(counts_all_mc)
counts_all_mc<-rbind(counts_all_mc,sc_annot)

#Filter to show only UT cells
counts_all_mc<-counts_all_mc[counts_all_mc$condition=="UT",]

#Create plot
g.2<-ggplot(counts_all_mc,aes(x=method,y=counts,fill=method))+
  geom_boxplot()+
  ylab("Number of (meta)cells per sample")+
  xlab("Method")+
  scale_y_log10()+
  scale_fill_discrete("Method")+
  theme(legend.position = "none",
        axis.title = element_text(size=14),
        axis.text = element_text(size=13),
        legend.title = element_text(size=13),
        legend.text= element_text(size=13))

# ------------------------------------------------------------------------------
# BLUEPRINT comparison
# ------------------------------------------------------------------------------

res<-read.table("metacell_general/eval_allmethods/perCondition_eval.tsv",header=TRUE)
res2<-read.table("metacell_general/eval_allmethods/singleCell_eval.tsv",header=TRUE)
res3<-read.table("metacell_general/eval_allmethods/sc_leiden_SCT_eval.tsv",header=TRUE)
#Parse method
res$method<-ifelse(grepl("leiden",res$File),"leiden",
                   ifelse(grepl("MetaCellar",res$File),"MetaCellaR","metacell"))
res2$method<-"singleCell"
res3$method<-ifelse(grepl("leiden",res3$File),"leiden_SCT","singleCell_SCT")

res<-rbind(res,res2,res3)
rm(res2,res3)

#Parse cutoff
res$cutoff<-as.numeric(stringi::stri_match(res$File,regex="cutoff(.*?)(_|\\.)")[,2])
res$cutoff<-paste0(as.character(res$cutoff*10),"%")

#Filter it to show only UT results and allGenes
res<-res[res$Condition=="UT" & res$Test == "allGenes",]

#Show only SCT results (also used later and noMetaCellaR)
res<-res[res$method %in% c("leiden_SCT","metacell","singleCell_SCT"),]
rename_methods<-setNames(c("Leiden","MetaCell","SingleCell"),
                         c("leiden_SCT","metacell","singleCell_SCT"))
res$method<-rename_methods[res$method]

g.3<-ggplot(res,aes(x=cutoff,y=Corr_corr,fill=method))+
  geom_bar(stat="identity",position="dodge")+
  ylab("Correlation with BLUEPRINT")+
  xlab("Genes stratified by x% expression in single cell")+
  scale_fill_discrete("Method")+
  theme(legend.position = "bottom",
        axis.title = element_text(size=14),
        axis.text = element_text(size=13),
        legend.title = element_text(size=13),
        legend.text= element_text(size=13))

g.bottom<-ggarrange(g.2,g.3,ncol=2,widths=c(0.4,0.6),
                    common.legend = TRUE,legend="bottom",
                    labels=c("b)","c)"))
g<-ggarrange(g.1,g.bottom,ncol=1,
             labels=c("a)",""))
ggsave(g,file="metacell_general/plots/metacell_overview_suppfigure.pdf",
       width=9,height=9)
ggsave(g,file="metacell_general/plots/metacell_overview_suppfigure.png",
       width=9,height=9)

