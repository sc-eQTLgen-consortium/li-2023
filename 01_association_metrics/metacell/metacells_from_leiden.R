# ------------------------------------------------------------------------------
# Implement own method to generate metacells based on leiden clustering
# Run leiden clustering separatley for each donor (run on Oelen v3, Monocytes)
# and use group cells that are part of the same cluster
# ------------------------------------------------------------------------------

library(Seurat)

#Load complete seurat object
seurat<-readRDS("../../seurat_objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds")
DefaultAssay(seurat)<-"SCT" #3000 most variable genes already identified

#Filter for monocytes
seurat<-seurat[,seurat$cell_type_lowerres=="monocyte"]

#Resolution for leiden clusters
leidenRes<-100
print(paste("Leiden resolution:",leidenRes))

type<-"SCT" #choose RNA or SCT
print(paste("Normalization:",type))

#Files with overall annotation and metacell matrix
annot_mc_all<-NULL
annot_mc_major_all<-NULL
metacellBulk_all<-NULL

#Iterate over all samples
samples<-levels(seurat$assignment)
for(donor in samples){
  
  print(paste("Processing donor:",donor))
  
  #Filter for the donor
  seurat_donor<-seurat[,seurat$assignment==donor]
  
  #Calculate PCA
  seurat_donor<-RunPCA(seurat_donor, verbose=FALSE)
  
  #Generate kNN graph and leidern clustering
  seurat_donor <- FindNeighbors(seurat_donor, dims = 1:20)
  seurat_donor <- FindClusters(seurat_donor, resolution = leidenRes, 
                               algorithm = 4, #4=Leiden
                               group.singletons=FALSE)
                          #don't assign all singletons to the nearest cluster 
  
  #Save metacell - cell annotation
  annot_mc<-data.frame(cluster=Idents(seurat_donor),
                       metacell=paste0("mc_",Idents(seurat_donor),"_",donor),
                       sample=donor,
                       cell=names(Idents(seurat_donor)),
                       row.names=NULL)
  annot_mc_all<-rbind(annot_mc_all,annot_mc)
  
  #Create pseudobulk
  #all(colnames(seurat_donor)==annot_mc$cell)
  if(type=="RNA"){
    metacellBulk <- t(apply(as.matrix(seurat_donor@assays$RNA@counts), 1, tapply, 
                            as.factor(annot_mc$cluster),
                            mean, na.rm=T))
  } else if (type=="SCT"){
    metacellBulk <- t(apply(as.matrix(seurat_donor@assays$SCT@counts), 1, tapply, 
                            as.factor(annot_mc$cluster),
                            mean, na.rm=T))    
  } else {
    stop(paste("Matrix type",type,"not known! Only RNA or SCT!"))
  }

  
  colnames(metacellBulk)<-paste0("mc_",1:ncol(metacellBulk),"_",donor)
  metacellBulk_all<-cbind(metacellBulk_all,metacellBulk)
  
  #Get majority annotation
  meta.data<-seurat_donor@meta.data
  meta.data$cell<-rownames(meta.data)
  meta.data<-merge(meta.data,annot_mc,
                   by.x="cell",by.y="cell")
  
  # Annotate each meta-cell to the most frequent condition
  timepoint.mc<-sapply(colnames(metacellBulk),
                       function(id) names(which.max(table(
                         meta.data$timepoint[
                           meta.data$metacell==id]))))
  
  #Save majority annotation
  annot_mc_major<-data.frame(metacell=names(timepoint.mc),
                             condition=unlist(timepoint.mc),
                             sample=donor,
                             row.names=NULL)
  
  annot_mc_major_all<-rbind(annot_mc_major_all,annot_mc_major)  
  
}



if(type=="RNA"){
  #Save per cell annotation
  write.table(annot_mc_all,file="annotations_metacell_leiden_perCell.tsv",sep="\t")
  write.table(annot_mc_major_all,file="annotations_mc_leiden_tp.tsv",sep="\t")
  
  #Save peudobulk counts
  saveRDS(metacellBulk_all, file="metacell_leiden.RDS")
} else if(type=="SCT"){
  write.table(annot_mc_all,file=paste0("annotations_metacell_leiden_SCT_perCell_",
                                       leidenRes,".tsv"),
              sep="\t")
  write.table(annot_mc_major_all,file=paste0("annotations_mc_leiden_SCT_tp_",
                                             leidenRes,".tsv"),
              sep="\t")
  
  #Save peudobulk counts
  saveRDS(metacellBulk_all, file=paste0("metacell_leiden_SCT_",
                                        leidenRes,".RDS"))
}


