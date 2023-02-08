require(Seurat)
require(slingshot)
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)

all <- readRDS('1M_v3_mediumQC_sct_celltyped_minimized_rnascaled.rds')
degenes <- read.table('degenes_monocyteUTX3hCA.txt')$V1
mono1 <- subset(x = all, subset = cell_type == 'mono 1')
mono1Ca <- subset(mono1, subset = (timepoint == 'UT') | (timepoint == 'X3hCA') | (timepoint == 'X24hCA'))
library(Matrix)
writeMM(GetAssayData(mono1Ca, assay='SCT', slot='data'), 
        "mono1Ca_allgenes.mtx")
write.table(as.matrix(mono1Ca[[]]), 'mono1Ca_allgenes.meta.csv', sep=",")
write.table(as.matrix(rownames(mono1Ca)), 'mono1Ca_allgenes.genes.txt')
mono1Ca_de3h <- subset(mono1, subset = (timepoint == 'UT') | (timepoint == 'X3hCA') | (timepoint == 'X24hCA'),
                  features = degenes)
# also select DE genes

# plot
pdf("pca_umap_sling_mono1CA_degenesUTX3h.pdf")
mono1Ca_de3h <- RunPCA(mono1Ca_de3h, npcs=10)
mono1Ca_de3h <- FindNeighbors(mono1Ca_de3h, verbose = FALSE, dims = 1:10)
mono1Ca_de3h <- FindClusters(mono1Ca_de3h, pc=1:10, algorithm = 2, random.seed = 256, resolution = 0.8)
mono1Ca_de3h <- RunUMAP(mono1Ca_de3h, dims = 1:10, reduction = "pca")

DimPlot(mono1Ca_de3h, reduction = "pca",
        group.by = "timepoint", pt.size = 0.5, label = TRUE, repel = TRUE)
ElbowPlot(mono1Ca_de3h, ndims=10)
DimPlot(mono1Ca_de3h, reduction = 'umap',
        group.by = "lane", pt.size = 0.5, label = TRUE, repel = TRUE)
DimPlot(mono1Ca_de3h, pt.size = 0.5, reduction = "umap", 
        group.by = "timepoint", label = TRUE)
DimPlot(mono1Ca_de3h, pt.size = 0.5, reduction = "umap", 
        group.by = "SCT_snn_res.0.8", label = TRUE)

# slingshot
mono1sling <- slingshot(Embeddings(mono1Ca_de3h, "umap"), clusterLabels = mono1Ca_de3h$SCT_snn_res.0.8, 
                       start.clus = 0, stretch = 0)
saveRDS(mono1Ca_de3h, 'mono1Ca_degenes.Rda')
saveRDS(mono1sling, 'mono1sling_degenes.Rda')
# load the expression data
mono1Ca_degenes <- readRDS('mono1Ca_degenes.Rda')
# load the slingshot
mono1sling <- readRDS('mono1sling_degenes.Rda')
pdf("evaluateK_chooseKnots.pdf")
mono1ca_matrix <- as.matrix(GetAssayData(mono1Ca_de3h, slot='data'))
icMat <- evaluateK(counts = mono1ca_matrix, 
                   sds = mono1sling, k = 3:10, 
                   nGenes = 200, verbose = T)
pdf("slingshot_pseudotime.pdf")
pseudotime <- slingPseudotime(mono1sling)

nc <- 2
nms <- colnames(pseudotime)
nr <- ceiling(length(nms)/nc)
par(mfrow = c(nr, nc))
for (i in nms) {
  ggplot(data.frame(pseudotime), aes(x=i, color=timepoint)) +
    geom_histogram(fill="white", alpha=0.5, position="identity")
}

ggplot_frame = data.frame(pseudotime)
ggplot_frame$timepoint <- mono1Ca_degenes[[]]$timepoint
ggplot(ggplot_frame, aes(x=curve1, color=timepoint)) +
  geom_histogram(fill="white", alpha=0.5, position="identity")
ggplot(ggplot_frame, aes(x=curve2, color=timepoint)) +
  geom_histogram(fill="white", alpha=0.5, position="identity")
ggplot(ggplot_frame, aes(x=curve3, color=timepoint)) +
  geom_histogram(fill="white", alpha=0.5, position="identity")
ggplot(ggplot_frame, aes(x=curve4, color=timepoint)) +
  geom_histogram(fill="white", alpha=0.5, position="identity")
dev.off()

#cellWeights <- slingCurveWeights(mono1sling)
#sce <- fitGAM(counts = GetAssayData(mono1Ca, slot='data'), 
#              pseudotime = pseudotime, cellWeights = cellWeights,
library(viridis)
pdf('slingshot_cells_in_different_linearges.pdf')
nc <- 2
nms <- colnames(pseudotime)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pseudotime[,i], breaks = 100)]
  plot(reducedDim(mono1sling), col = colors, pch = 16, cex = 0.5, main = i)
  lines(mono1sling, lwd = 2, col = 'black', type = 'lineages')
}

dev.off()


library(SCORPIUS)
mono1Ca_degenes <- readRDS('mono1Ca_degenes.Rda')
#pdf('SCORPIUS_plots.pdf')
expression <- t(as.matrix(GetAssayData(mono1Ca_degenes, slot='data')))
group_name <- factor(as.character(mono1Ca_degenes[[]]$timepoint))
# try with PCA
#pdf('scorpius_pca.pdf')
#pearson_space <- reduce_dimensionality(expression, "pearson")
#pearson_traj <- infer_trajectory(pearson_space)
#draw_trajectory_plot(pearson_space, group_name, pearson_traj$path, contour = TRUE)
#dev.off()

space <- reduce_dimensionality(expression, "spearman")
traj <- infer_trajectory(space)
saveRDS(space, 'scorpius_space.rds')
saveRDS(traj, 'scorpius_traj.rds')
# save traj#time in tsv
write.table(as.matrix(traj$time), file='scorpius_trajtime.tsv', sep = '\t')
write.table(as.matrix(traj$path), file='scorpius_trajpath.tsv', sep = '\t')

# load scorpius results
space <- readRDS('scorpius_space.rds')
traj <- readRDS('scorpius_traj.rds')
draw_trajectory_plot(space, group_name, traj$path, contour = TRUE)
histogram_data <- data.frame("time" = matrix(unlist(traj$time), nrow=length(traj$time), byrow=T),
                             row.names=names(traj$time))
histogram_data$timepoint <- mono1Ca_degenes[[]]$timepoint

ggplot(histogram_data, aes(x=time, color=timepoint)) +
  geom_histogram(fill="white", alpha=0.5, position="identity")

# draw_trajectory_heatmap(space, traj$time, progression_group=group_name)
pdf('scorpius_heatmap.pdf')
gimp <- gene_importances(
  expression, 
  traj$time, 
  num_permutations = 10, 
  num_threads = 8, 
  ntree = 10000,
  ntree_perm = 1000
) 
saveRDS(gimp, 'scorpius_gimp.rds')
gimp$qvalue <- p.adjust(gimp$pvalue, "BH", length(gimp$pvalue))
gene_sel <- gimp$gene[gimp$qvalue < .05]
expr_sel <- scale_quantile(expression[,gene_sel])

modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = T) # needs more RAM than 50G
draw_trajectory_heatmap(expr_sel, traj$time, group_name, modules)
dev.off()

