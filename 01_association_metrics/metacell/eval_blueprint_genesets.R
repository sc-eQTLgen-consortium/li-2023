# ------------------------------------------------------------------------------
# Compare correlation from BLUEPRINT with correlation from metacells/single cell
# for different expression tresholds
# Here: shown for leiden metacells, calculation done the same way for original
#       metacells and single cell
# ------------------------------------------------------------------------------

library(reticulate) # to read the single cell data (numpy)
library(data.table)
library(ggplot2) # only required if plotting=TRUE

np <- import("numpy")

#Iterate over the list of different gene sets
corr_files<-c("metacell_general/leiden_metacells/correlation_r_leiden_SCT_cutoff08.tsv",
              "metacell_general/leiden_metacells/correlation_r_leiden_SCT_cutoff06.tsv",
              "metacell_general/leiden_metacells/correlation_r_leiden_SCT_cutoff04.tsv",
              "metacell_general/leiden_metacells/correlation_r_leiden_SCT_cutoff02.tsv")

res_file<-"metacell_general/eval_allmethods/sc_leiden_SCT_eval.tsv"

plotting<-TRUE

#Write column headers
write.table(data.frame("Condition","Num_pairs","Corr_corr","Test","File"),
            file=res_file,quote=FALSE,sep="\t", 
            row.names=FALSE, col.names = FALSE)

#Load the large blueprint data set
path<-"blueprint/allGenePairs_BlueprintScMonocytes_GeneGeneCorrelationComparison.pairwiseSpearman."
corr.blue.vals <- np$load(paste0(path,"npy"))
corr.blue.vals<-corr.blue.vals[,1]

corr.blue<-fread(paste0(path,"genePairs.txt"),header=FALSE)

#Split into gene1 and gene2
corr.blue$gene1<-sapply(corr.blue$V1,function(s) strsplit(s,"/")[[1]][1])
corr.blue$gene2<-sapply(corr.blue$V1,function(s) strsplit(s,"/")[[1]][2])
corr.blue$V1<-NULL
corr.blue$corr.blue<-corr.blue.vals
rm(corr.blue.vals)

for(cfile in corr_files){
  
  print(paste("Processing:",cfile))
  
  #Load corr.mc.r
  corr.mc.r<-read.table(cfile)
  
  #Filter Blueprint matrix
  corr.genes<-unique(c(corr.mc.r$Gene1,corr.mc.r$Gene2))
  corr.blue.subset<-corr.blue[gene1 %in% corr.genes & 
                         gene2 %in% corr.genes]
  
  #Order correctly so that gene1 smaller than gene2
  corr.blue.subset$swap<-ifelse(corr.blue.subset$gene1 < corr.blue.subset$gene2, 
                                corr.blue.subset$gene1,corr.blue.subset$gene2)
  corr.blue.subset$gene2<-ifelse(corr.blue.subset$gene1 < corr.blue.subset$gene2, 
                                 corr.blue.subset$gene2,corr.blue.subset$gene1)
  corr.blue.subset$gene1<-corr.blue.subset$swap
  corr.blue.subset$swap<-NULL
  colnames(corr.blue.subset)<-c("Gene1","Gene2","Correlation.blue")
  
  #Merge everything
  corr.mc.r<-merge(corr.mc.r,corr.blue.subset,by=c("Gene1","Gene2"))
  corr.mc.r<-reshape2::melt(corr.mc.r,id.vars=c("Gene1","Gene2","Correlation.blue"))
  colnames(corr.mc.r)[4:5]<-c("Condition","Correlation")
  
  #Correlation for all genes
  corr.corr<-sapply(unique(corr.mc.r$Condition), function(tp)  
    cor(corr.mc.r$Correlation.blue[corr.mc.r$Condition==tp],
        corr.mc.r$Correlation[corr.mc.r$Condition==tp],
        method="pearson",use = "pairwise.complete.obs"))
  
  res<-data.frame(condition=unique(corr.mc.r$Condition),
                  num.pairs=as.vector(table(corr.mc.r$Condition)),
                  corr.corr,
                  test="allGenes",
                  file=cfile)
  
  write.table(res, file=res_file,quote=FALSE,sep="\t",
              append=TRUE,row.names=FALSE, col.names = FALSE)

  corr.mc.r<-corr.mc.r[! is.na(corr.mc.r$Correlation.blue),]
  if(plotting){
    corr.mc.r$class<-ifelse(corr.mc.r$Correlation.blue>0,
                            "positive","negative")
    
    g<-ggplot(corr.mc.r,aes(x=Correlation.blue,y=Correlation,color=class))+
      geom_point()+facet_wrap(~Condition,ncol=3)+
      xlab("Correlation Blueprint")+ylab("Correlation MC")
    ggsave(g,file=paste0("metacell_general/eval_allmethods/plots/",
                         "comp_corr_blue_",strsplit(cfile,"/")[[1]][2],".png"))
  }
  
  #Correlation for genes with positive correlation
  corr.mc.r.pos<-corr.mc.r[corr.mc.r$Correlation.blue>0,]
  corr.corr<-sapply(unique(corr.mc.r.pos$Condition), function(tp)  
    cor(corr.mc.r.pos$Correlation.blue[corr.mc.r.pos$Condition==tp],
        corr.mc.r.pos$Correlation[corr.mc.r.pos$Condition==tp],
        method="pearson",use = "pairwise.complete.obs"))
  
  res<-data.frame(condition=unique(corr.mc.r.pos$Condition),
                  num.pairs=as.vector(table(corr.mc.r.pos$Condition)),
                  corr.corr,
                  test="posGenes",
                  file=cfile)
  
  write.table(res, file=res_file,quote=FALSE,sep="\t",
              append=TRUE,row.names=FALSE, col.names = FALSE)
  
  #Correlation for genes with negative correlation
  corr.mc.r.neg<-corr.mc.r[corr.mc.r$Correlation.blue<0,]
  corr.corr<-sapply(unique(corr.mc.r.neg$Condition), function(tp)  
    cor(corr.mc.r.neg$Correlation.blue[corr.mc.r.neg$Condition==tp],
        corr.mc.r.neg$Correlation[corr.mc.r.neg$Condition==tp],
        method="pearson",use = "pairwise.complete.obs"))
  
  res<-data.frame(condition=unique(corr.mc.r.neg$Condition),
                  num.pairs=as.vector(table(corr.mc.r.neg$Condition)),
                  corr.corr,
                  test="negGenes",
                  file=cfile)
  
  write.table(res, file=res_file,quote=FALSE,sep="\t",
              append=TRUE,row.names=FALSE, col.names = FALSE)
}
