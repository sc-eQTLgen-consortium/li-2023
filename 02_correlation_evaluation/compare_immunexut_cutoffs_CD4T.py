# ---------------------------------------------------------------------------------------
# Compare correlation between ImmuNexUT and single cell (Oelen v3 dataset)
# for different thresholds (number of cells expressing the gene),
# implemented for UT and CD4T cells here
# Input: seurat objects with Oelen v3 dataset and normalized ImmuNexUT counts
# Output: csv file with the correlation between ImmuNexUT and Oelen v3 for each threshold
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

# load ImmuNexuT object
counts = pd.read_csv('imd_paper_rna_data/norm_count/Naive_CD4_norm_count.txt',sep="\t")
immunexut_genes = counts.index.values
counts = counts.transpose()

# method to select genes above a certain nonzero ratio
def select_gene_nonzeroratio(df, ratio):
    nonzerocounts = np.count_nonzero(df.values, axis=0)/df.shape[0]
    selected_genes = df.columns[nonzerocounts>ratio]
    return selected_genes

# generate a set of thresholds that should be tested (start with stricter thresholds)
thresholds = [i/10 for i in range(1,10)]
thresholds.reverse()

f_out = open("co-expression_indivs_combined/immunexut_cutoff_eval_CD4T.txt", "w")
f_out.write("threshold,ngenes,corr_pearson\n")

# iterate over all thresholds
for th in thresholds:
    
    #select all genes within the threshold
    selected_genes = select_gene_nonzeroratio(celltype_data, th)
    
    # filter genes that are not in ImmuNexUT
    selected_genes = list(set(selected_genes) & set(immunexut_genes))
    
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
    
    #c alculate correlation ImmuNexUT
    input_df_ImmuNexUT = counts[selected_genes]
    input_data = spearmanr(input_df_ImmuNexUT, axis=0)[0]
    input_data_uppertria = input_data[np.triu_indices_from(input_data, 1)]

    corrs_df_ImmuNexUT = pd.DataFrame({'BULK': input_data_uppertria},
                    index=gene_pairs)
    
    # sorting both is not necessary here
    #all(corrs_df.index == corrs_df_ImmuNexUT.index)
    
    # calculate correlation between datasets and save results
    corr_data = pearsonr(corrs_df.UT, corrs_df_ImmuNexUT.BULK)[0]
    
    # save results
    f_out.write(f"{th},{len(selected_genes)},{corr_data}\n")

# close file
f_out.close()

