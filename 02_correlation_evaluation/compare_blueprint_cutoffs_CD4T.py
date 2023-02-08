# ---------------------------------------------------------------------------------------
# Compare correlation between Blueprint and single cell (Oelen v3 dataset)
# for different thresholds (number of cells expressing the gene),
# implemented for UT and CD4T cells here
# Input: seurat objects with Oelen v3 dataset and precalculated Blueprint correlation
#        for all possible gene pairs
# Output: csv file with the correlation between Blueprint and Oelen v3 for each threshold
# ---------------------------------------------------------------------------------------

from scipy.stats import spearmanr, pearsonr
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
from time import time
import os
import re

# load scanpy object (Oelen v3 dataset)
alldata = sc.read_h5ad('seurat_objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.SCT.h5ad')

# filter for CD4+ T cells and UT cells
alldata = alldata[alldata.obs.cell_type_lowerres=='CD4T']
alldata = alldata[alldata.obs.timepoint=='UT'].copy() #copy to not create only a view object

celltype_data = pd.DataFrame(data=alldata.X.toarray(),
                                 index=alldata.obs.index,
                                 columns=alldata.var.index)

# load Blueprint object
bp_corr = np.load('blueprint_data/tcel_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.spearmanR.npy',mmap_mode="r")

bp_corr_genes = []
f= open('blueprint_data/tcel_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.spearmanR.cols.txt','r')
for line in f.readlines():
    bp_corr_genes.append(line.rstrip())

# method to select genes above a certain nonzero ratio
def select_gene_nonzeroratio(df, ratio):
    nonzerocounts = np.count_nonzero(df.values, axis=0)/df.shape[0]
    selected_genes = df.columns[nonzerocounts>ratio]
    return selected_genes

# generate a set of thresholds that should be tested (start with stricter thresholds)
thresholds = [i/10 for i in range(1,10)]
thresholds.reverse()

f_out = open("co-expression_indivs_combined/blueprint_cutoff_eval_CD4T.txt", "w")
f_out.write("threshold,ngenes,corr_pearson\n")

# iterate over all thresholds
for th in thresholds:
    
    #select all genes within the threshold
    selected_genes = select_gene_nonzeroratio(celltype_data, th)
    
    # filter genes that are not in Blueprint
    selected_genes = list(set(selected_genes) & set(bp_corr_genes))
    
    print(f"Number of selected genes for {th}: {len(selected_genes)}")
    
    gene_pairs = []
    for i,gene1 in enumerate(selected_genes):
        for j in range(i+1, len(selected_genes)):
            if gene1 < selected_genes[j]:
                gene_pairs.append(';'.join([gene1, selected_genes[j]]))
            else:
                gene_pairs.append(';'.join([selected_genes[j],gene1]))

    # calculate correlation single cell
    input_df = celltype_data[selected_genes]
    input_data = spearmanr(input_df, axis=0)[0]
    input_data_uppertria = input_data[np.triu_indices_from(input_data, 1)]

    corrs_df = pd.DataFrame({'UT': input_data_uppertria},
                            index=gene_pairs)
    
    # filter blueprint and order it the same way as the single cell object
    filter_bp_genes = [gene in selected_genes for gene in bp_corr_genes]
    bp_corr_filtered = bp_corr[filter_bp_genes][:,filter_bp_genes]
    bp_uppertria = bp_corr_filtered[np.triu_indices_from(bp_corr_filtered, 1)]
    
    # get genes from the blueprint object
    bp_corr_genes_filtered = [gene for gene in bp_corr_genes if gene in selected_genes]
    gene_pairs_bp=[]
    for i,gene1 in enumerate(bp_corr_genes_filtered):
        for j in range(i+1, len(bp_corr_genes_filtered)):
            if gene1 < bp_corr_genes_filtered[j]:
                gene_pairs_bp.append(';'.join([gene1, bp_corr_genes_filtered[j]]))
            else:
                gene_pairs_bp.append(';'.join([bp_corr_genes_filtered[j],gene1]))
                
    corrs_df_bp = pd.DataFrame({'BP': bp_uppertria},
                            index=gene_pairs_bp)
    
    # sort both and combine them
    corrs_df = corrs_df.sort_index()
    corrs_df_bp = corrs_df_bp.sort_index()
    #all(corrs_df.index == corrs_df_bp.index)
    
    # calculate correlation between datasets and save results
    corr_data = pearsonr(corrs_df.UT, corrs_df_bp.BP)[0]
    
    # save results
    f_out.write(f"{th},{len(selected_genes)},{corr_data}\n")

# close file
f.close()
