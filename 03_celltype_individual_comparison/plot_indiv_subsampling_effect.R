# ------------------------------------------------------------------------------
# Plot effect of number of cells on correlation between individuals
# Input: pairwise comparison of all individuals (Pearson correlation) per cell type
#        and for different numbers of cells 
#        (subsampling and calculation done in correlation_subsampling.py)
# Output: violin plot showing trend
# ------------------------------------------------------------------------------

library(ggplot2)

theme_set(theme_bw())

suffix<-"v3"
#suffix<-"v2"

color_coding <- list()
color_coding[["CD4+ T"]] <- "#2E9D33"
color_coding[["CD8+ T"]] <- "#126725"
color_coding[["Monocyte"]] <- "#EDBA1B"
color_coding[["NK"]] <- "#E64B50"
color_coding[["B"]] <- "#009DDB"
color_coding[["DC"]] <- "#965EC8"

#Full cell type names as reported in the paper
cell_types_corrected<-setNames(c("CD4+ T","CD8+ T","Monocyte","NK","DC","B"),
                               c("CD4T","CD8T","monocyte","NK","DC","B"))

#Load results
res<-read.csv(paste0("co-expression_indivs_subsampled/",
                     "correlation_individuals_subsampled_1M_",suffix,".csv"))
res$X<-NULL

res$celltype<-cell_types_corrected[res$celltype]

#Filter out some values to make it more visible
res<-res[res$cell_num %in% seq(25,500,50),]

res$cell_num<-as.factor(res$cell_num)

g<-ggplot(res,aes(x=cell_num,y=corr,fill=celltype))+
  geom_violin(position = position_dodge(0.9)) +
  xlab("Subsampled number of cells per individual")+
  ylab("Correlation between individuals")+
  ylim(0,1)+
  scale_fill_manual("Cell type",values=color_coding)+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        legend.title = element_text(size=13),
        legend.text = element_text(size=13),
        legend.position=c(0.9,0.2))+
  guides(fill=guide_legend(nrow=3,byrow=FALSE))
print(g)
ggsave(g,file=paste0("co-expression_indivs_subsampled/plots/subsampling_1M_",
                     suffix,"_filtered.pdf"),
       width=14,height=4)
