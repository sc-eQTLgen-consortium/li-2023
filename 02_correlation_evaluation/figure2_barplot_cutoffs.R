# ------------------------------------------------------------------------------
# Create barplot of correlation dependency on expression cutoff
# for Oelen v3 and ImmuNexUT / Blueprint
# (only plotting in R, calculation done with python script)
# Input: correlation comparison between Blueprint and Oelen v3 dataset 
#        (precalculated in compare_blueprint_cutoffs_CD4T.py) and between
#        ImmuNexUT and Oelen v3 dataset (precalculated in 
#        compare_immunexut_cutoffs_CD4T.py)
# Output: two barplots, one for Blueprint comparsion and one for 
#         ImmuNexUT comparison
# ------------------------------------------------------------------------------

library(ggplot2)
library(RColorBrewer)

theme_set(theme_bw())

################################################################################
# Plot for ImmuNexUT (main Figure 2c)
################################################################################

vals<-read.table("co-expression_indivs_combined/immunexut_cutoff_eval_CD4T.txt",
                 sep=",",header=TRUE)

vals$threshold<-as.factor(vals$threshold)
g<-ggplot(vals,aes(x=threshold,y=corr_pearson,fill=ngenes))+
  geom_bar(stat="identity")+
  geom_text(aes(x = threshold, y = corr_pearson / 2, label = ngenes,
            color=ifelse(ngenes<1000,'white','black')),size=5)+
  scale_color_manual(values=c("black","white"))+
  xlab("Expression cutoff")+
  ylab("Correlation between Oelen (v3)\nand ImmuNexUT")+ylim(0,1)+
  scale_fill_distiller("Number of\ngenes",palette="YlOrBr")+
  theme(legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        legend.title=element_text(size=13),
        legend.text=element_text(size=12))+
  guides(color=FALSE)
print(g)
ggsave(g,file="co-expression_indivs_combined/plots/eval_immunexut_cutoff.pdf",
       width=6,height=5)


################################################################################
# Plot for Blueprint (Supplement)
################################################################################

vals<-read.table("co-expression_indivs_combined/blueprint_cutoff_eval_CD4T.txt",
                 sep=",",header=TRUE)

vals$threshold<-as.factor(vals$threshold)
g<-ggplot(vals,aes(x=threshold,y=corr_pearson,fill=ngenes))+
  geom_bar(stat="identity")+
  geom_text(aes(x = threshold, y = corr_pearson / 2, label = ngenes,
                color=ifelse(ngenes<1000,'white','black')),size=5)+
  scale_color_manual(values=c("black","white"))+
  xlab("Expression cutoff")+
  ylab("Correlation between Oelen (v3)\nand BLUEPRINT")+ylim(0,1)+
  scale_fill_distiller("Number of\ngenes",palette="YlOrBr")+
  theme(legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        legend.title=element_text(size=13),
        legend.text=element_text(size=12))+
  guides(color=FALSE)
print(g)
ggsave(g,file="co-expression_indivs_combined/plots/eval_blueprint_cutoff.png",
       width=6,height=5)

