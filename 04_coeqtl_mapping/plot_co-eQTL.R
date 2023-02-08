############################################################################################################################
# Code Author: Dylan de Vries
# Name: plot_co-eQTL.R
# Function: Plot co-eQTLs
############################################################################################################################
#
# Libraries
#
############################################################################################################################
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(gridExtra)

############################################################################################################################
#
# Functions
#
############################################################################################################################
# Name: get.expression.data
# Function: Get the expression data and calculate the co-expression for plotting purposes
# Input:
#   Name 	            Type          Description
#   sample				character	  sample name
#	cell.type			character	  cell type to get the data for
# 	gene1 				character 	  first gene to get data for
# 	gene2 				character 	  second gene to get data for
# 	genotype			character 	  the genotype of the co-eQTL for this sample 
#
# Output:
# 	A list with two data frames of one sample. The first is for making the boxplots and the second for the personalized expression regression plot
get.expression.data <- function(sample, cell.type, gene1, gene2, genotype){
	sample.gene1.expression <- data@assays$SCT@data[gene1, rownames(data@meta.data[data@meta.data$cell_type_lowerres == cell.type & data@meta.data$assignment == sample,])]
	sample.gene2.expression <- data@assays$SCT@data[gene2, rownames(data@meta.data[data@meta.data$cell_type_lowerres == cell.type & data@meta.data$assignment == sample,])]

	sample.co.expression <- cor(sample.gene1.expression, sample.gene2.expression, method="spearman")
	expr.plot.data <- data.frame(gene1.expression=sample.gene1.expression, gene2.expression=sample.gene2.expression, sample=sample, genotype=genotype)
	plot.data <- list(sample.co.expression, expr.plot.data)
	return(plot.data)
}

# Name: prepare.plot.data
# Function: Combine the data of all samples into data.frames
# Input:
#   Name 	            Type          Description
# 	gene1 				character 	  first gene to get data for
# 	gene2 				character 	  second gene to get data for
# 	SNP.name			character 	  the rs-ID for the co-eQTL SNP
#	cell.type			character	  cell type to get the data for
#
# Output:
# 	A list with two data frames. The first is for making the boxplots and the second for the personalized expression regression plot
prepare.plot.data <- function(gene1, gene2, SNP.name, cell.type){
	co.expressions <- c()
	genotypes <- c()
	expr.plot.data <- data.frame(gene1.expression=numeric(0), gene2.expression=numeric(0), sample=character(0), genotype=character(0))
	for (sample in samples){
		genotypes <- c(genotypes, genotypes_all[SNP.name, sample])
		plot.data <- get.expression.data(sample, cell.type, gene1, gene2, genotypes_all[SNP.name, sample])
		expr.plot.data <- rbind(expr.plot.data, plot.data[[2]])
		co.expressions <- c(co.expressions, plot.data[[1]])
	}
	plot.data <- data.frame(co.expression=co.expressions, sample=samples, genotype=genotypes)
	combined.plot.data <- list(plot.data, expr.plot.data)
	return(combined.plot.data)
}

# Name: plot.co.eQTL.boxplot
# Function: Make a plot for the co-eQTL
# Input:
#   Name 	            Type          Description
#   plot.data 			data.frame 	  the data for the boxplot
#   expr.plot.data		data.frame 	  the data for the expression regression plot
# 	gene1 				character 	  first gene to get data for
# 	gene2 				character 	  second gene to get data for
# 	SNP.name			character 	  the rs-ID for the co-eQTL SNP
#	cell.type			character	  cell type to get the data for
#	meta.z				numeric		  meta z-score
#	QTL.type			character	  indicates whether it's amongst the strongest, middle or weakest co-eQTLs
#	QTL.type.index		character	  the index of the co-eQTL within its type
#
# Output:
# 	A list with two data frames. The first is for making the boxplots and the second for the personalized expression regression plot
plot.co.eQTL.boxplot <- function(plot.data, expr.plot.data, gene1, gene2, SNP.name, cell.type, meta.z, QTL.type, QTL.type.index){
	genotype.colors <- c("#57a350", "#fd7600", "#383bfe", "white")
	names(genotype.colors) <- c("0/0", "0/1", "1/1", "white")

	sample.color <- c(colorRampPalette(c("#9efc95", "#57a350"))(length(which(plot.data$genotype=="0/0"))),
		colorRampPalette(c("#fabb84", "#fd7600"))(length(which(plot.data$genotype=="0/1"))),
		colorRampPalette(c("#acadfc", "#383bfe"))(length(which(plot.data$genotype=="1/1"))))
	names(sample.color) <- c(as.character(plot.data$sample[plot.data$genotype == "0/0"]), as.character(plot.data$sample[plot.data$genotype == "0/1"]), as.character(plot.data$sample[plot.data$genotype == "1/1"]))

	expr.plot <- ggplot(expr.plot.data, aes(x=gene1.expression, y=gene2.expression, fill=sample, color=sample)) + geom_point(size=0.5) +
		geom_smooth(method = "lm", fullrange = T, se=F) +
		scale_fill_manual(values=sample.color) +
		scale_color_manual(values=sample.color) +
		xlab(paste0(gene1, " expression")) +
		ylab(paste0(gene2, " expression")) +
		ggtitle(paste0(SNP.name, " effect on ", gene1, " - ", gene2, "\nco-expression")) +
		guides(fill=FALSE, color=FALSE) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7), panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))

	box.plot <- ggplot(plot.data) + geom_boxplot(aes(x=genotype, y=co.expression, fill=genotype), outlier.shape=NA, alpha=0.6) + 
		geom_quasirandom(aes(x=genotype, y=co.expression, color=genotype, fill="white"), pch=21, size=2, alpha=1, dodge.width=0.4, alpha=0.6) +
		scale_fill_manual(values=genotype.colors) + 
		scale_color_manual(values=genotype.colors) +
		xlab("Genotype") +
		ylab(paste0(gene1, " - ", gene2, " co-expression")) +
		ggtitle(paste0(SNP.name, " co-eQTL\n", QTL.type, " ", QTL.type.index)) +
		guides(fill=FALSE, color=FALSE) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7), panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))

	pdf(paste0("/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/co-eQTLs/plots/", cell.type, "/", cell.type, "_co-eQTL_", SNP.name, "_", gene1, "-", gene2, ".pdf"))
	grid.arrange(expr.plot, box.plot, ncol=2)
	dev.off()
}

############################################################################################################################
#
# Main code
#
############################################################################################################################
data <- readRDS("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/seurat_objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.rds")
vcf <- fread('/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/genotypes/LL_trityper_plink_converted.vcf.gz')
target.QTLs <- read.table("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/output/filtered_results/UT_monocyte/coeqtls_fullresults.sig.tsv.gz", header=T, sep="\t", stringsAsFactors=F)
target.QTLs <- target.QTLs[order(abs(target.QTLs$MetaPZ), decreasing=T),]
genotypes_all <- as.data.frame(vcf[, 10:ncol(vcf)])
rownames(genotypes_all) <- vcf$ID

#Get the 10 strongest, 10 middling and 10 weakest of the input co-eQTLs
QTL.selection <- target.QTLs[c(1:10, floor(nrow(target.QTLs)/2):(floor(nrow(target.QTLs)/2)+10), (nrow(target.QTLs)-10):nrow(target.QTLs)),]
samples <- unique(data@meta.data$assignment)

for (QTL.index in 1:nrow(QTL.selection)){
	print(QTL.index)
	if (QTL.index <= 10){
		type <- "strong"
		QTL.type.index <- QTL.index
	} else if (QTL.index <= 20){
		type <- "medium"
		QTL.type.index <- QTL.index - 10
	} else {
		type <- "poor"
		QTL.type.index <- QTL.index - 20
	}
	genes <- unlist(strsplit(QTL.selection$Gene[QTL.index], ";"))
	combined.plot.data <- prepare.plot.data(genes[1], genes[2], QTL.selection$SNP[QTL.index], "monocyte")
	plot.data <- combined.plot.data[[1]]
	expr.plot.data <- combined.plot.data[[2]]

	plot.co.eQTL.boxplot(plot.data, expr.plot.data, genes[1], genes[2], QTL.selection$SNP[QTL.index], "monocyte", QTL.selection$MetaPZ[QTL.index], type, QTL.type.index)
}
