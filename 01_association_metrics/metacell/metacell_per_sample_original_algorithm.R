# ------------------------------------------------------------------------------
# Metacell algorithm (original) run for each sample separately
# with Oelen v3 dataset (Monocytes)
# Take all 200 variable genes 
# (but removing for each sample the ones with too low coverage)
#
# Remarks: 
# * metacells uses a "data base", the processed files are not loaded
#   directly in the workspace
# * in the function mcell_mc_from_coclust_balanced the parameters
#   K and min_mc_size can be used to change the "size" of meta cells,
#   but too small is not recommended
#
# ------------------------------------------------------------------------------

library(metacell)
library(SingleCellExperiment)
library(ggplot2)
library(optparse)

#Parse arguments
option_list = list(
  make_option(c("-K","--coclustK"), default="20",
              help="Parameter K of mcell_mc_from_coclust_balanced",
              type="integer"),
  make_option(c("-m","--coclustMinMcSize"), default="10",
              help="Parameter min_mc_size of mcell_mc_from_coclust_balanced",
              type="integer")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#Parameters to change the granularity of the samples
mc_coclustK<-opt$coclustK
mc_coclustMin_size<-opt$coclustMinMcSize

print(paste("Running meta cells with following parameters (changing granularity):",
            "K",mc_coclustK,"min_mc_size",mc_coclustMin_size))

#Create directories to save the meta cells and the correlation results
mainDir<-paste0("metacells_K",mc_coclustK,"_minCells",mc_coclustMin_size)
if(!dir.exists(mainDir)) 
  dir.create(mainDir)

setwd(mainDir)

##########################################################################
# Important note: a lot of the parameter are manged over tgconfig
# see: https://github.com/tanaylab/tgconfig

# To check all set parameters
#tgconfig::get_package_params('metacell')

#Set number of cores (otherwise I get issues for very small data sets)
tgconfig::set_param('mc_cores', 8, 'metacell')

#Issues with downsampling matrix => set parameter for downsampling lower 
#(probably problem as matrix is too sparse)
#tgconfig::set_param("scm_n_downsamp_gstat",300,'metacell')

########################################################################

#Option to create additional plots for visualization
allPlots<-FALSE
fileType<-"seurat"

#Load the h5ad object and convert it to a single cell object
if(fileType=="h5ad"){
  
  library(reticulate) # to load h5ad object
  library(zellkonverter) # to convert h5ad object to single cell object
  
  sc<-import("scanpy")
  adata<-sc$read("../../seurat_objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.SCT.h5ad")
  #Filter for monocytes
  adata<-adata[adata$obs$cell_type_lowerres=="monocyte"]$copy()
  sce<-AnnData2SCE(adata)
  rm(adata)
  
  #Convert assay name to counts (before called X)
  assayNames(sce)<-c("counts")

  #Alternatively read a seurat object  
} else if (fileType=="seurat"){
  
  library(Seurat)
  
  seurat<-readRDS("../../seurat_objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds")
  
  #Filter for monocytes
  seurat<-seurat[,seurat$cell_type_lowerres=="monocyte"]
  
  #Convert into a single cell object
  sce <- as.SingleCellExperiment(seurat,assay="RNA")
  rm(seurat)
  
} else {
  stop("File type not known!")
}

#Select all remaining samples
samples<-as.character(unique(sce@colData$assignment))

#Save the results
metacell.allsamples<-NULL
annotations.allsamples<-NULL
annotations.percell<-NULL
#Calculate meta-cell algorithm for each sample separatly
#Also option to calculate correlation, but currently not done
for(sample in samples){
  
  print(paste("Processing sample:",sample))
  sce.sample.full<-sce[,sce@colData$assignment == sample]

  #Create data base directory
  if(!dir.exists("database")){
    dir.create("database/")
  } else {
    do.call(file.remove, list(list.files("database/", full.names = TRUE)))
  }
  scdb_init("database/", force_reinit=T)
  
  #Upload SCE object
  #Filter for genes with at least 4 counts (preliminary before the real filtering downstream 
  #to reduce the calculation burden)
  mat<-scm_import_sce_to_mat(sce.sample.full[rowSums(counts(sce.sample.full))>3,])
  scdb_add_mat(sample, mat)
  
  #Create a directory for figures
  if(!dir.exists("figs")) dir.create("figs/")
  scfigs_init("figs/")
  
  #Create a gset for generating the knn graph 
  mcell_add_gene_stat(gstat_id="stat", mat_id=sample)
  mcell_gset_filter_varmean(gset_id="sample_feat", gstat_id="stat", T_vm=0.08, force_new=T)
  #Sampled coverage of at least T_tot and threshold for the third highest UMI count > T_top3
  mcell_gset_filter_cov(gset_id = "sample_feat", gstat_id="stat", T_tot=100, T_top3=2)
  #Check generated gene set
  gset<-scdb_gset("sample_feat")
  print(paste("Number of selected genes:",length(gset@gene_set)))
  
  #Create the knn graph based on correlation
  mcell_add_cgraph_from_mat_bknn(mat_id=sample,
                                 gset_id = "sample_feat",
                                 graph_id="sample_graph",
                                 K=50,
                                 dsamp=T)
  
  #Resample cells from the graph to robustly define groups
  mcell_coclust_from_graph_resamp(
    coc_id="sample_coc500",
    graph_id="sample_graph",
    min_mc_size=20,
    p_resamp=0.75,
    n_resamp=500)
  
  #Remark the size of the meta cells can be influenced by the paramters 
  #K and min_mc_size (for both is true: the smaller, the more cells ...)  
  mcell_mc_from_coclust_balanced(
    coc_id="sample_coc500",
    mat_id= sample,
    mc_id= paste0(sample,"_mc"),
    K=mc_coclustK,
    min_mc_size=mc_coclustMin_size, 
    alpha=2)
  
  #Plotting outlier (only possible for small groups)
  if(allPlots){
    mcell_plot_outlier_heatmap(mc_id=paste0(sample,"_mc"), 
                               mat_id = sample, T_lfc=3) 
  }

  #Split and filter metacells using dbscan and outlier gene detection
  mcell_mc_split_filt(new_mc_id=paste0(sample,"_mc_f"),
                      mc_id=paste0(sample,"_mc"),
                      mat_id=sample,
                      T_lfc=3, plot_mats=F)
  
  ##Selecting marker genes automatically
  mcell_gset_from_mc_markers(gset_id="sample_markers", mc_id=paste0(sample,"_mc_f"))
  mc_colorize_default(paste0(sample,"_mc_f"))
  
  #Creating a heatmap of genes and metacells 
  #(also not really well visible with too many cells)
  if(allPlots){
    mcell_mc_plot_marks(mc_id=paste0(sample,"_mc_f"), gset_id="sample_markers", 
                        mat_id=sample)
  }

  #Create graph layout
  mcell_mc2d_force_knn(mc2d_id=paste0(sample,"_2dproj"),
                       mc_id=paste0(sample,"_mc_f"), graph_id="sample_graph")
  #Plotting also not really interesting for two large 
  mcell_mc2d_plot(mc2d_id=paste0(sample,"_2dproj"))
  
  #Save it again as a h5ad object to compare the results
  #So far no direct exporting function found, therefore processing the object myself
  #See https://tanaylab.github.io/metacell/reference/tgMCCov-class.html
  sce_meta<-scdb_mc(paste0(sample,"_mc_f"))
  
  #Meta cell annotations
  mc.annot<-data.frame(metaCell=sce_meta@mc)
  mc.annot$cell<-rownames(mc.annot)
  rownames(mc.annot)<-NULL
  mc.annot<-rbind(mc.annot,
                  data.frame(metaCell=0,
                             cell=sce_meta@outliers))
  
  annotations.percell<-rbind(annotations.percell,
                             mc.annot)
    
  #Check distributions between cell types and stimulation results
  annotations<-sce.sample.full@colData
  #annotations$cell<-rownames(annotations)
  annotations<-merge(mc.annot,annotations,by.x="cell",by.y="bare_barcode_lane")
    
  perMetacell<-as.data.frame(table(annotations$metaCell))
  
  #Plot only timepoint for now
  freqs<-as.data.frame(table(annotations$metaCell,
                             annotations$timepoint))
  freqs<-merge(freqs,perMetacell,by="Var1",suffixes=c(".spc",".mc"))
  freqs$Fraction<-freqs$Freq.spc/freqs$Freq.mc

  g<-ggplot(freqs,aes(x=as.factor(Var1),y=Fraction,fill=Var2))+
      geom_bar(stat="identity")+
      xlab("Meta cell (0=Outlier)")+
      scale_fill_discrete(name = "Time point")+
      ggtitle(paste("Cell number in total:",sum(freqs$Freq.spc)))
  ggsave(g,filename=paste0("figs/barplot_time_ct_",sample,".png"))
  
  #Create a pseudobulk object with the meta-cell annotation (without outliers)
  mc.annot<-mc.annot[mc.annot$metaCell>0,]
  sc.counts<-counts(sce.sample.full)[,mc.annot$cell]
  all(colnames(sc.counts)==mc.annot$cell)
  
  mc.annot$metaCell<-as.factor(paste0(sample,"_mc_",mc.annot$metaCell))
  mc.pseudobulk<- t(apply(sc.counts, 1, tapply, mc.annot$metaCell,
               sum, na.rm=T))
  
  #Normalize to 10,000 per metacell
  libSize<-colSums(mc.pseudobulk)
  mc.pseudobulk<-t(t(mc.pseudobulk)/libSize*10000)
  
  metacell.allsamples<-cbind(metacell.allsamples,mc.pseudobulk)
  
  #Create a majority annotation for each metacell
  annotations$metaCell<-paste0(sample,"_mc_",annotations$metaCell)
  timepoint.mc<-sapply(colnames(mc.pseudobulk),
                       function(id) names(which.max(table(
    annotations$timepoint[
      annotations$metaCell==id]))))
  
  annot.mc<-data.frame(metacell=names(timepoint.mc),
                       timepoint=timepoint.mc,
                       sample=sample)
  
  #Add how many cells where part of the meta-cell
  perMetacell$Var1<-paste0(sample,"_mc_",perMetacell$Var1)
  colnames(perMetacell)<-c("metacell","cell.count")
  annot.mc<-merge(annot.mc,perMetacell,by="metacell")
  annotations.allsamples<-rbind(annotations.allsamples,annot.mc)
  
}

#Save results
write.table(annotations.allsamples,file="annotations_metacell.tsv",sep="\t")
write.table(annotations.percell,file="annotations_singlecell_metacell.tsv",sep="\t")
saveRDS(metacell.allsamples,file="pseudobulk_metacell.RDS")

#Delete the database directory
unlink("database",recursive = TRUE)

