##############################################################################
# Calculate correlation for each cell type for the van Blokland v3 dataset,
# timpeoint 6-8 weeks after admission, merging all individuals
##############################################################################

#from scipy.stats import t, norm
from scipy.stats import spearmanr
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
from time import time
import os
import re

# specify if the gene selection was done before and is passed in a file
gene_selection_file = False

# set working directory (to shorten path length)
os.chdir('./')

# load scanpy object
prefix_results = Path('co-expression_indivs_combined/stemi/version3')
# test stemi v2
alldata = sc.read_h5ad('seurat_objects/cardio.integrated.20210301.stemiv3.h5ad')

def select_gene_nonzeroratio(df, ratio):
    nonzerocounts = np.count_nonzero(df.values, axis=0)/df.shape[0]
    selected_genes = df.columns[nonzerocounts>ratio]
    return selected_genes

# extract timepoint from timepoint - stimulation annotation
def get_time(x):
    if x == 'UT':
        return x
    else:
        pattern = re.compile(r'\d+h')
        return re.findall(pattern, x)[0]


# extract timepoint from timepoint - stimulation annotation
observations = alldata.obs.copy()
observations['timepoint_id_celltype'] = [f'{item[0]}_{item[1]}' for item
                                           in observations[['timepoint.final', 'cell_type_lowerres']].values]

celltypes = ['B', 'CD4T', 'CD8T', 'monocyte', 'DC', 'NK']
for celltype in celltypes:
    if not os.path.isdir(prefix_results/celltype):
        os.mkdir(prefix_results/celltype)
    starttime = time()
    print(celltype)
    specific = alldata[alldata.obs.cell_type_lowerres==celltype]
    celltype_data = pd.DataFrame(data=specific.X.toarray(),
                                 index=specific.obs.index,
                                 columns=specific.var.index)

    # get the set of gene pairs
    specific_obs = observations[observations['cell_type_lowerres']==celltype]
    for condition in observations['timepoint.final'].unique():
        # filter for the condition
        celltype_condition_data = celltype_data[specific_obs['timepoint.final']==condition]
        
        # take either tsv file with selected genes or filter genes after a nonzero rate
        if gene_selection_file:
            selected_genes = pd.read_csv(' /genelists_tp_union/expressed_gene_'+celltype+'.tsv')
            selected_genes =selected_genes["genes"].tolist()
        else:
            selected_genes = select_gene_nonzeroratio(celltype_condition_data, 0.5)

        print(f"Number of selected genes for {celltype} {condition}: {len(selected_genes)}")

        gene_pairs = []
        for i,gene1 in enumerate(selected_genes):
            for j in range(i+1, len(selected_genes)):
                gene_pairs.append(';'.join([gene1, selected_genes[j]]))

        input_df = celltype_condition_data[selected_genes]
        input_data = spearmanr(input_df, axis=0)[0]
        input_data_uppertria = input_data[np.triu_indices_from(input_data, 1)]

        corrs_df = pd.DataFrame(data=input_data_uppertria,
                                columns=[f'{condition}'],
                                index=gene_pairs)

        corrs_df.to_csv(prefix_results/celltype/f'{celltype}_{condition}_correlation.csv')

        #Filter for 0.3 correlation cutoff
        corrs_df = corrs_df[corrs_df[condition]>0.3]
        corrs_df.to_csv(prefix_results/celltype/f'{celltype}_{condition}_correlation_03filtered.csv')

        print(f"Finished {celltype} with time {time() - starttime}")
