# ------------------------------------------------------------------------------
# Check Pearson correlation between individuals (per cell type)
# and combine it with correlation levels in the cell type 
# (as both are below each other in the final figure)
# Input: 
# 1) correlation matrices per cell type generated with 
#    correlation_timepoint_combined_indivs.py
#    (for correlation levels in each cell type)
# 2) correlation matrices per individual and cell type 
#    (for comparison of individuals)
# Output: plot and summary as output text
#  -----------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(dplyr)

theme_set(theme_bw())

path<-"coeqtl_mapping/input/individual_networks/UT/"

cell_types<-c("B","CD4T","CD8T","DC","monocyte","NK")

#Full cell type names as reported in the paper
cell_types_corrected<-setNames(c("CD4+ T","CD8+ T","Monocyte","NK","DC","B"),
                               c("CD4T","CD8T","monocyte","NK","DC","B"))

#Get standard color scale
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#Evaluate both Oelen v2 and v3 dataset
for(dataset in c("onemillionv2","onemillionv3")){
  
  #Evaluate correlation distribution in the cell type
  all_corrs<-NULL
  tp<-"UT"
  for(ct in cell_types){
    
    #Load correlation values
    if(dataset=="onemillionv3"){
      corr_ct<-fread(paste0("co-expression_indivs_combined/",ct,"/",
                            ct,"_",tp,"_correlation.csv"))
    } else {
      corr_ct<-fread(paste0("co-expression_indivs_combined/one_million_version2/",
                            ct,"/",ct,"_",tp,"_correlation.csv"))     
    }
    colnames(corr_ct)[2]<-"corr"
    
    #Get absolute correlation
    corr_ct$corr<-abs(corr_ct$corr)
    all_corrs<-rbind(all_corrs,
                     data.frame(level=c("<0.05",">0.05",">0.1",">0.2",">0.3"),
                                values=c(sum(abs(corr_ct$corr<0.05)),
                                         sum(abs(corr_ct$corr)>0.05 &
                                               abs(corr_ct$corr)<0.1),
                                         sum(abs(corr_ct$corr)>0.1 &
                                               abs(corr_ct$corr)<0.2),
                                         sum(abs(corr_ct$corr)>0.2 & 
                                               abs(corr_ct$corr)<0.3),
                                         sum(abs(corr_ct$corr)>0.3)),
                                freq=c(mean(abs(corr_ct$corr<0.05)),
                                       mean(abs(corr_ct$corr)>0.05 &
                                              abs(corr_ct$corr)<0.1),
                                       mean(abs(corr_ct$corr)>0.1 &
                                              abs(corr_ct$corr)<0.2),
                                       mean(abs(corr_ct$corr)>0.2 & 
                                              abs(corr_ct$corr)<0.3),
                                       mean(abs(corr_ct$corr)>0.3)),
                                ct,tp))
  }
  
  #Check general distribution of "highly correlated genes"
  high_corr<-all_corrs[all_corrs$level %in% c(">0.3",">0.2",">0.1"),]
  high_corr<-high_corr%>%group_by(ct)%>%
    summarize(high_freq=sum(freq))%>%
    as.data.frame()
  median(high_corr$high_freq)  
  
  #Get comparison between cell types
  summary<-NULL
  for(ct in cell_types){
    #Load file with all correlation values per individual and cell type
    #each individual one column
    if(dataset=="onemillionv3"){
      corr<-fread(paste0(path,dataset,"/UT_",ct,".genesnonzero0.5.coefs.gz"))
    } else {
      corr<-fread(paste0(path,dataset,"/UT_",ct,".genesnonzero0.5.coefs.tsv.gz"))
    }
  
    #Get Pearson correlation between individual-specific correlations
    gene_pairs<-corr$V1
    corr$V1<-NULL
    indiv_corr<-cor(corr,method="pearson")
    
    #Melt the upper triangle
    tmp<-reshape2::melt(indiv_corr)
    tmp$Var1<-as.character(tmp$Var1)
    tmp$Var2<-as.character(tmp$Var2)
    tmp<-tmp[tmp$Var1 < tmp$Var2,]
    
    tmp$ct<-ct
    
    summary<-rbind(summary,tmp)
  }
  
  #Replace cell type names
  all_corrs$ct<-cell_types_corrected[all_corrs$ct]
  summary$ct<-cell_types_corrected[summary$ct]
  
  #Sort cell types according to their highly correlated genes
  sorting<-all_corrs[all_corrs$level==">0.3",]
  
  #Sort cell type colors
  colors_cts<-gg_color_hue(6)
  colors_cts<-colors_cts[order(sorting$freq)]
  
  sorting<-sorting[order(sorting$freq),]
  
  #Barplot showing the general correlation distribution in the cell type
  all_corrs$level<-factor(all_corrs$level,levels=c("<0.05",">0.05",">0.1",">0.2",">0.3"))
  all_corrs$ct<-factor(all_corrs$ct,levels=sorting$ct)
  colors_bars<-brewer.pal(n = 6, "YlGnBu")
  g.1<-ggplot(all_corrs,aes(x=ct,y=freq,fill=level))+
    geom_bar(stat="identity")+
    xlab("Cell type")+ylab("Fraction correlated genes")+
    scale_fill_manual("Absolute\ncorrelation",values=colors_bars[2:6])+
    theme(legend.position="bottom",
          axis.title.y = element_text(size=13.5),
          axis.title.x = element_blank(),
          axis.text = element_text(size=12),
          legend.title=element_text(size=11),
          legend.text=element_text(size=10.5))
  
  #Violin plot showing differences between individuals in the cell type
  summary$ct<-factor(summary$ct,levels=sorting$ct)
  g.2<-ggplot(summary,aes(x=ct,fill=ct,y=value))+
    geom_violin()+
    geom_boxplot(width = 0.15, outlier.shape = NA)+
    ylim(0,1)+
    xlab("Cell type")+
    ylab("Correlation between individuals")+
    scale_fill_manual(values=colors_cts)+
    theme(legend.position = "none",
          axis.title.y = element_text(size=13.5),
          axis.title.x = element_blank(),
          axis.text = element_text(size=12))

    
  g<-ggarrange(g.1,g.2,ncol=1,align="hv")
  ggsave(g,file=paste0("co-expression_indivs_combined/plots/corr_ct_indiv_",
                       dataset,".pdf"),
         width=5,height=6)
  
  #Get median correlation in each cell type
  med_corr<-summary%>%
    group_by(ct)%>%
    summarise(median(value,na.rm=TRUE))
  
  print(dataset)
  print(med_corr)
  
}
