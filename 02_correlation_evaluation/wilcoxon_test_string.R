# ------------------------------------------------------------------------------
# Compare if correlated pairs from single cell (Oelen v3, CD4+ T cells)
# and bulk (ImmuNexUT, naive CD4+ T cells) are enriched in STRING database
# (Using the same strategy as in CRISPR validation with 
# Wilcoxon Rank Sum Test)
# ------------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(RColorBrewer)

theme_set(theme_bw())

#Get colors
cols_brewer <- c("grey78",brewer.pal(n = 3, "Set2")[2])

cond<-"UT"

plot_list<-NULL
for(data_type in c("sc","ImmuNexUT")){
  
  print(data_type)
  
  # Load single cell correlation (cell type specific)
  if(data_type == "sc"){
    
    ct<-"CD4T" #alternative "CD8T"
    
    corr_ct<-fread(paste0("co-expression_indivs_combined/",ct,"/",ct,"_",cond,
                          "_correlation.csv"))
    corr_ct$gene1<-gsub(";.*","",corr_ct$V1)
    corr_ct$gene2<-gsub(".*;","",corr_ct$V1)
    
    corr_ct$swap<-ifelse(corr_ct$gene1 > corr_ct$gene2,corr_ct$gene1,corr_ct$gene2)
    corr_ct$gene1<-ifelse(corr_ct$gene1 > corr_ct$gene2,corr_ct$gene2,corr_ct$gene1)
    corr_ct$gene2<-corr_ct$swap
    corr_ct$swap<-NULL
  } else if (data_type == "ImmuNexUT"){
    
    ct<-"Naive_CD4"
    
    #Read ImmuNexUT data (already preprocessed correctly)
    corr_ct<-fread(paste0("imd_paper_rna_data/correlation/",
                          ct,"_correlation.txt"))
    corr_ct$V1<-NULL
    colnames(corr_ct)[3]<-"UT"
  }
  
  expressed_genes<-union(corr_ct$gene1,corr_ct$gene2)
  
  corr_ct$UT<-abs(corr_ct$UT)
  
  #Read STRING data base
  string<-fread("additional_files/STRING-network.csv")
  string<-string[string$Gene1 %in% expressed_genes &
                   string$Gene2 %in% expressed_genes,]  
    
  string$swap<-ifelse(string$Gene1 > string$Gene2,string$Gene1,string$Gene2)
  string$Gene1<-ifelse(string$Gene1 > string$Gene2,string$Gene2,string$Gene1)
  string$Gene2<-string$swap
  string$swap<-NULL
  
  #Combine with correlation values
  string$is_string<-TRUE
  corr_ct<-merge(corr_ct,string,by.x=c("gene1","gene2"),
                 by.y=c("Gene1","Gene2"),all.x=TRUE)
  corr_ct$is_string[is.na(corr_ct$is_string)]<-FALSE
  
  wt<-wilcox.test(corr_ct$UT[corr_ct$is_string],corr_ct$UT[!corr_ct$is_string],
                  paired=FALSE,alternative="greater")
  print(wt$p.value)
  
  if(data_type=="sc"){
    data_type<-"Oelen (v3)"
  }
  
  g<-ggplot(corr_ct,aes(x=is_string,y=UT, fill=is_string))+
    geom_violin()+
    geom_boxplot(width = 0.15, outlier.shape = NA)+
    xlab("Gene pair in STRING network")+
    ylab("Absolute correlation")+
    ylim(c(0,1))+
    ggtitle(paste(data_type,"dataset"))+
    annotate("text",x=1.5,y=0.9,label=paste0("p =",
            format(wt$p.value,digits=2)),size=4.5)+
    scale_fill_manual(values=cols_brewer)+
    theme(legend.position = "none",
          plot.title=element_text(size=15),
          axis.title=element_text(size=14),
          axis.text=element_text(size=12))
  plot_list<-c(plot_list,list(g))
}

g<-ggarrange(plotlist=plot_list,ncol=2,labels=c("a","b"))
ggsave(g,file=paste0("compare_with_string/plots/string_wilcoxon_combined.pdf"),
       width=7,height=4)