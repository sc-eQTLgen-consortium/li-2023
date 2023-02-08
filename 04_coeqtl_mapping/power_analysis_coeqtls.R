################################################################################
# Evaluate how number of triplets decreases the power to detect
# co-expression QTLs by calculating the power dependent on the sample size
# (N=173), heritability (Rsq: 0.1-0.3), Bonferroni multiple testing correction,
# and different number of tests
# The number of tests is estimated based on different expression cutoffs for
# the Oelen v3 dataset, assuming that all pairwise combinations are tested for
# all genes above the respective cutoff and one SNP per pair
# Input: Seurat object with data from Oelen v3
# Output: line plot visualizing power for different number of tests
################################################################################

library(Seurat)
library(scPower)
library(ggplot2)

theme_set(theme_bw())

################################################################################
# Getting expression distribution for Oelen v3 dataset
################################################################################

#Load complete seurat object
seurat<-readRDS("seurat_objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds")

#Filter for monocytes and UT timepoint
seurat<-seurat[,seurat$cell_type_lowerres == "monocyte"]
seurat<-seurat[,seurat$timepoint == "UT"]

#Calculate for each gene the non-zero ratio
nonzero_ratio<-rowMeans(as.matrix(seurat@assays$SCT@counts)>0)

#Get cumulative ratio
thresholds<-seq(0,1,0.05)
num_genes<-sapply(thresholds,function(i)sum(nonzero_ratio>i))

#Save results in a file
nonzero_count<-data.frame(nonzero_ratio=thresholds,num_genes)

################################################################################
# Performing power calculation
################################################################################

#Samples in meta-analysis
nSamples<-173

bonfLevel<-function(nTests){
  return(0.05/nTests)
}

#Test different heritabilities
Rsq<-seq(0.1,0.3,0.05)

#Number tests
nonzero_count$genepairs<-nonzero_count$num_genes*(nonzero_count$num_genes-1)/2

res<-NULL
for(her in Rsq){
  for(i in 1:(nrow(nonzero_count)-1)){
    res<-rbind(res,
               data.frame(her,
                          numTests=nonzero_count$genepairs[i],
                          cutoff=nonzero_count$nonzero_ratio[i],
                          power=scPower:::power.eqtl.ftest(her,
                                        bonfLevel(nonzero_count$genepairs[i]),
                                                           nSamples)))
  }
}

#Plot results
g<-ggplot(res,aes(x=numTests,y=power,color=as.factor(her)))+
  geom_line()+
  scale_color_discrete("Heritability")+
  scale_x_log10()+
  xlab("Number tests")+ylab("Power")
print(g)
ggsave(g,file="power_calculation/power_effect_nonzeroratio.png",
       height=5,width=6)

