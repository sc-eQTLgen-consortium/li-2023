#############################################################################################
# Calculate correlation for each pairwise gene combiation in each cell type and the UT timepoint,
# where both genes are expressed in at least 50% of the cells
# merging all individuals for Oelen v2 and v3 dataset (specified in parameter version2)
# Input: seurat objects with Oelen v2 and v3 dataset
# Output: correlation values as csv files (one per cell type)
############################################################################################

#from scipy.stats import t, norm
from scipy.stats import spearmanr
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
from time import time
import os
import re

# specify if Oelen v2 (version2=True) or v3 (version2=False) dataset is used
version2 = True

# set result path
if version2:
    prefix_results = Path('co-expression_indivs_combined/one_million_version2/')
else:
    prefix_results = Path('co-expression_indivs_combined/')

# load scanpy object
if version2:
    alldata = sc.read_h5ad('seurat_objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.sct.h5ad')
else:
    alldata = sc.read_h5ad('seurat_objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.SCT.h5ad')
    
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
observations['time_merged'] = [get_time(item) for item in observations['timepoint']]
observations['timepoint_id_celltype'] = [f'{item[0]}_{item[1]}' for item
                                           in observations[['time_merged', 'cell_type_lowerres']].values]

# iterate over each cell type
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

    # select only UT cells
    for condition in ['UT']:

        # filter for the condition
        celltype_condition_data = celltype_data[specific_obs.time_merged==condition]
        
        #filter genes after a nonzero rate of at least 0.5
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


        print(f"Finished {celltype} with time {time() - starttime}")
