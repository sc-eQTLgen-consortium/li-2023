##############################################################################
# Calculate correlation for each cell type for the van Blokland v2 dataset,
# timpeoint 6-8 weeks after admission, merging all individuals
##############################################################################

from scipy.stats import spearmanr
import numpy as np
import pandas as pd
from pathlib import Path
from time import time
import os


def select_gene_nonzeroratio(df, ratio):
    nonzerocounts = np.count_nonzero(df.values, axis=0)/df.shape[0]
    selected_genes = df.columns[nonzerocounts>ratio]
    return selected_genes

# set working directory (to shorten path length)
os.chdir('./')

# load scanpy object
prefix_results = Path('co-expression_indivs_combined/stemi/')
stemi_data = pd.read_csv('seurat_objects/stemi_v2_monocyte.csv.gz', compression='gzip', sep=' ', index_col=0).T
stemi_meta = pd.read_csv('seurat_objects/stemi_v2_monocyte.meta.csv', sep=' ', index_col=0)

for condition in stemi_meta['timepoint.final'].unique():
    starttime = time()
    # filter for the condition
    celltype_condition_data = stemi_data[stemi_meta['timepoint.final']==condition]
    # take either tsv file with selected genes or filter genes after a nonzero rate
    selected_genes = select_gene_nonzeroratio(celltype_condition_data, 0.5)
    print(f"Number of selected genes for stemi {condition}: {len(selected_genes)}")
    # get gene pair names
    gene_pairs = []
    for i,gene1 in enumerate(selected_genes):
        for j in range(i+1, len(selected_genes)):
            gene_pairs.append(';'.join([gene1, selected_genes[j]]))
    # get gene-gene correlations
    input_df = celltype_condition_data[selected_genes]
    input_data = spearmanr(input_df, axis=0)[0]
    input_data_uppertria = input_data[np.triu_indices_from(input_data, 1)]
    corrs_df = pd.DataFrame(data=input_data_uppertria,
                            columns=[f'{condition}'],
                            index=gene_pairs)
    corrs_df.to_csv(prefix_results/f'monocyte_{condition}_correlation.csv')
    #Filter for 0.3 correlation cutoff
    corrs_df = corrs_df[corrs_df[condition]>0.3]
    corrs_df.to_csv(prefix_results/f'monocyte_{condition}_correlation_03filtered.csv')
    print(f"Finished {condition} with time {time() - starttime}")
