# ------------------------------------------------------------------------------
# Gene selections for metacell evaluation: generate files for different
# gene subsets (expressed in x% until (x+20)% of the cells with x between 20-80)
# for Oelen v3 dataset, Monocytes
# This allows threshold dependent evaluation for BLUEPRINT comparison. 
# See details downstream scripts metacell_general_correlation_tp.R, 
# single_cell_correlation_tp.R, eval_blueprint_genesets.R
# ------------------------------------------------------------------------------

library(Seurat)

#Load complete seurat object
seurat<-readRDS("seurat_objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds")
DefaultAssay(seurat)<-"RNA"

#Filter for monocytes
seurat<-seurat[,seurat$cell_type_lowerres=="monocyte"]

#Selected cutoffs
cutoffs<-c(1,0.8,0.6,0.4,0.2)

#Split into lists dependent on expression cutoff
exprGenes.singleCell<-rowSums(as.matrix(seurat@assays$RNA@counts)>0)/ncol(seurat)

print(paste("Number of genes expressed in at least 50% of cells:",
            sum(exprGenes.singleCell>=0.5)))

for(i in 1:(length(cutoffs)-1)){
  gene.subset<-rownames(seurat)[exprGenes.singleCell<=cutoffs[i] &
                                  exprGenes.singleCell>cutoffs[i+1]]
  print(paste("Number of genes with expression cutoff",
              cutoffs[i+1],":",length(gene.subset)))
  
  write.table(gene.subset,
              file=paste0("metacell_general/eval_allmethods/gene_lists/",
                          "mono_expr_genes_cut_",cutoffs[i+1],".txt"),
              row.names = FALSE,col.names = FALSE,quote=FALSE)
}

