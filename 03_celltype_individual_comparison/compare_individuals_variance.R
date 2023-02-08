# ------------------------------------------------------------------------------
# Combine gene pair variance across individuals 
# for Oelen (v2) and (v3) in one plot (taking Z scores)
# Input: correlation matrices per individual and cell type 
#        (for comparison of individuals)
# Output: plot and summary as output text
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

theme_set(theme_bw())

path<-"coeqtl_mapping/input/individual_networks/UT/"

cell_types<-c("B","CD4T","CD8T","DC","monocyte","NK")

#Full cell type names as reported in the paper
cell_types_corrected<-setNames(c("CD4+ T","CD8+ T","Monocyte","NK","DC","B"),
                               c("CD4T","CD8T","monocyte","NK","DC","B"))

g_list<-NULL
#Evaluate both Oelen v2 and v3 dataset
for(dataset in c("onemillionv2","onemillionv3")){
  corr_summary<-NULL
  for(ct in cell_types){
    
    #Correlation values
    if(dataset=="onemillionv2"){
      corr<-fread(paste0(path,dataset,"/UT_",ct,".genesnonzero0.5.zscores.tsv.gz"))
    } else {
      corr<-fread(paste0(path,dataset,"/UT_",ct,".genesnonzero0.5.zscores.gz"))
    }
    
    gene_pairs<-corr$V1
    corr$V1<-NULL
    
    #Set Inf values to NA to remove them afterwards
    corr<-as.matrix(corr)
    corr[is.infinite(corr)]<-NA
    
    #Get mean and variance for each gene pair (drop NA and Inf values from calculation)
    corr_summary<-rbind(corr_summary,
                   data.frame(ct,
                              gene_pairs,
                              var=apply(corr,1,var,na.rm=TRUE),
                              mean=apply(corr,1,mean,na.rm=TRUE)))
  }
  
  #Check frequency of "highly variable" genes
  print(paste("Frequency of highly variable genes for",dataset))
  tmp<-corr_summary%>%
    group_by(ct)%>%
    summarise(freq_high_var=mean(var>2,na.rm=TRUE))
  print(tmp)
  print(median(tmp$freq_high_var))
  
  #Replace cell type names
  corr_summary$ct<-cell_types_corrected[corr_summary$ct]
  
  g<-ggplot(corr_summary,aes(x=var,color=ct))+
    geom_density()+
    xlab("Correlation variance across individuals")+
    ylab("Density")+
    ylim(0,2)+
    xlim(0,15)+
    scale_color_discrete("Cell type")+
    theme(axis.title=element_text(size=10),
          axis.text=element_text(size=9),
          legend.title=element_text(size=10),
          legend.text=element_text(size=10))
  
  g_list<-c(g_list,list(g))
}

g<-ggarrange(plotlist=g_list,ncol=2,common.legend = TRUE,
             legend="bottom",labels=c("a","b"))
ggsave(g,file="co-expression_indivs_combined/plots/per_indivual_var_zscores_combined.pdf",
       width=8,height=4)