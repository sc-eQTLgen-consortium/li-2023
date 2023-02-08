# ------------------------------------------------------------------------------
# Calculate correlation per timepoint (and sample if stated) from the
# original single cell dataset (Oelen v3 dataset, Monocytes) for different
# gene sets (split dependent on gene expression cutoff) for comparison with
# metacells (see corresponding files create_genesets.R,
# metacell_general_correlation_tp.R and eval_blueprint_genesets.R)
# Input: Seurat object, file with selected genes
# Output: files with correlation values (r-values and p-values)
# ------------------------------------------------------------------------------

library(Seurat)
library(Hmisc) #for fast calculation of correlation
library(optparse)

#Parse arguments
option_list = list(
  make_option(c("-g","--selectedGenes"), 
              default="benchmark/celltypes/gene_expressed_over_hald_cells.txt",
              help="path to list with selected genes"),
  make_option(c("-s","--perSample"),action="store_true",
              default=FALSE,
              help="Shall the evaluation be done for each sample separatly"),
  make_option(c("-t","--type"),
              default="RNA",
              help="Use either RNA count matrix (RNA) or SCT count matrix (SCT)."),
  make_option(c("-o","--outputFile"),
              default="timepoint_monocytes",
              help="Suffix of the output files")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

pathSelectedGenes<-opt$selectedGenes
perSample<-opt$perSample
matrixType<-opt$type
outputSuffix<-opt$outputFile

print(paste("Evaluating single cell data for gene set:"))
print(pathSelectedGenes)

print(paste("Evaluating each sample individually:", perSample))

#Load complete seurat object
seurat<-readRDS("seurat_objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds")

#Filter for monocytes
seurat<-seurat[,seurat$cell_type_lowerres=="monocyte"]

#Select which genes shall be chosen for correlation (same as for single cell)
selected.genes<-read.table(pathSelectedGenes,
                           header=FALSE)
seurat<-seurat[selected.genes$V1,]

#Result data frame (correlation and pvalues)
corr.df<-NULL
pval.df<-NULL

correlationRes<-function(meta_counts,colName){
  
  #Be carefull: rcorr does not work with less than 5 samples
  corr.mc<-rcorr(t(meta_counts), type="spearman")
  #corr.mc<-cor(t(meta_counts), method="spearman")
  
  #Create a pairwise data frame for the correlation
  corr.pairs.mc<-as.data.frame(as.table(corr.mc$r),
                               stringsAsFactors = FALSE)
  corr.pairs.mc<-corr.pairs.mc[corr.pairs.mc$Var1<corr.pairs.mc$Var2,]
  colnames(corr.pairs.mc)<-c("Gene1","Gene2",colName)
  
  #Create a pairwise data frame for the pvalue
  corr.pairs.pval<-as.data.frame(as.table(corr.mc$P),
                                 stringsAsFactors = FALSE)
  corr.pairs.pval<-corr.pairs.pval[corr.pairs.pval$Var1<corr.pairs.pval$Var2,]
  colnames(corr.pairs.pval)<-c("Gene1","Gene2",colName)
  
  return(list(corr.pairs.mc,corr.pairs.pval))
}

if(perSample){
  
  #Calculate correlation for each sample
  samples<-unique(seurat$assignment)
  for(sample in samples){
    print(paste("Processing sample:",sample))

    for(timepoint in unique(seurat$timepoint[seurat$assignment==sample])){

      print(paste("Calculate correlation for timepoint",timepoint))

      meta_counts<-seurat[,seurat$assignment==sample &
                            seurat$timepoint==timepoint]
      if(matrixType=="RNA"){
        meta_counts<-as.matrix(meta_counts@assays$RNA@counts)
      } else if (matrixType=="SCT"){
        meta_counts<-as.matrix(meta_counts@assays$SCT@counts)
      } else {
        stop(paste("Matrix type",matrixType,"not known! Only RNA or SCT!"))
      }
      tmp<-correlationRes(meta_counts,colName = paste0(timepoint,"-",sample))
      corr.pairs.mc<-tmp[[1]]
      corr.pairs.pval<-tmp[[2]]

      #Concatinate the sample - timepoint pairs
      if(is.null(corr.df)){
        corr.df<-corr.pairs.mc
        pval.df<-corr.pairs.pval
      } else {
        corr.df<-merge(corr.df,corr.pairs.mc,by=c("Gene1","Gene2"),
                       all=TRUE)
        pval.df<-merge(pval.df,corr.pairs.pval,by=c("Gene1","Gene2"),
                       all=TRUE)
      }
    }
  }

} else {

  for(timepoint in unique(seurat$timepoint)){
    
    print(paste("Calculate correlation for timepoint",timepoint))
    
    meta_counts<-seurat[,seurat$timepoint==timepoint]
    if(matrixType=="RNA"){
      meta_counts<-as.matrix(meta_counts@assays$RNA@counts)
    } else if (matrixType=="SCT"){
      meta_counts<-as.matrix(meta_counts@assays$SCT@counts)
    } else {
      stop(paste("Matrix type",matrixType,"not known! Only RNA or SCT!"))
    }
    
    tmp<-correlationRes(meta_counts,colName = timepoint)
    corr.pairs.mc<-tmp[[1]]
    corr.pairs.pval<-tmp[[2]]
    
    #Concatinate the sample - timepoint pairs
    if(is.null(corr.df)){
      corr.df<-corr.pairs.mc
      pval.df<-corr.pairs.pval
    } else {
      corr.df<-merge(corr.df,corr.pairs.mc,by=c("Gene1","Gene2"),
                     all=TRUE)
      pval.df<-merge(pval.df,corr.pairs.pval,by=c("Gene1","Gene2"),
                     all=TRUE)
    }
    
  }
}

write.table(corr.df,
            file=paste0("metacell_general/eval_allmethods/corr_sc/correlation_r_",
                        outputSuffix,".tsv"),
            sep="\t",quote=FALSE)
write.table(pval.df,
            file=paste0("metacell_general/eval_allmethods/corr_sc/correlation_pval_",
                        outputSuffix,".tsv"),
            sep="\t",quote=FALSE)