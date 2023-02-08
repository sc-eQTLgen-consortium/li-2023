######################################################################
# Calculate correlation for each cell type for van der Wijst dataset,
# merging all individuals
######################################################################

from scipy.stats import spearmanr
import numpy as np
import pandas as pd
from pathlib import Path
from time import time
import scanpy as sc
import os


def select_gene_nonzeroratio(df, ratio):
    nonzerocounts = np.count_nonzero(df.values, axis=0)/df.shape[0]
    selected_genes = df.columns[nonzerocounts>ratio]
    return selected_genes


def select_gene_variances(df, ratio):
    variances = np.var(df.values, axis=0)/df.shape[0]
    var_thres = np.percentile(variances, ratio)
    # var_thres = np.nanmedian(variances)
    selected_genes = df.columns[variances>var_thres]
    print(selected_genes[:5])
    return selected_genes


def get_genename(df, mapping):
    df['genename'] = [mapping.get(geneid) for geneid in df.index]
    df = df.dropna(subset=['genename']).drop_duplicates(subset=['genename'])
    df = df.set_index('genename')
    return df


# set working directory (to shorten path length)
os.chdir('./')
gene_selection_file = False

# load scanpy object
prefix_results = Path('co-expression_indivs_combined/ng_updated_version')
# test stemi v2
# alldata = sc.read_h5ad('seurat_objects/pilot3_subsetted_celltypes_final_ensemble_converted_samples.h5ad')
alldata = sc.read_h5ad('seurat_objects/pilot3_seurat3_200420_sct_azimuth.h5ad')

# extract timepoint from timepoint - stimulation annotation
celltype_maping = {'CD4 T': 'CD4T', 'CD8 T': 'CD8T', 'Mono': 'monocyte', 'DC': 'DC', 'NK':'NK', 'other T': 'otherT', 'other': 'other', 'B':'B'}
alldata.obs['cell_type_mapped_to_onemillion'] = [celltype_maping.get(name) for name in alldata.obs['predicted.celltype.l1']]
observations = alldata.obs.copy()
celltypes = [item for item in observations['cell_type_mapped_to_onemillion'].unique() if not pd.isnull(item)]
print(celltypes)
for celltype in celltypes:
    if not os.path.isdir(prefix_results / celltype):
        os.mkdir(prefix_results / celltype)
    starttime = time()
    print(celltype)
    specific = alldata[alldata.obs['cell_type_mapped_to_onemillion'] == celltype]
    celltype_data = pd.DataFrame(data=specific.X.toarray(),
                                 index=specific.obs.index,
                                 columns=specific.var.index)
    print(celltype_data.shape)
    # get the set of gene pairs
    specific_obs = observations[observations['cell_type_mapped_to_onemillion'] == celltype]
    # filter for the condition
    celltype_condition_data = celltype_data
    # take either tsv file with selected genes or filter genes after a nonzero rate
    if gene_selection_file:
        selected_genes = pd.read_csv(' /genelists_tp_union/expressed_gene_' + celltype + '.tsv')
        selected_genes = selected_genes["genes"].tolist()
    else:
        selected_genes = select_gene_nonzeroratio(celltype_condition_data, 0.5)
    print(f"Number of selected genes for {celltype} : {len(selected_genes)}")
    gene_pairs = []
    for i, gene1 in enumerate(selected_genes):
        for j in range(i + 1, len(selected_genes)):
            gene_pairs.append(';'.join([gene1, selected_genes[j]]))
    input_df = celltype_condition_data[selected_genes]
    input_data = spearmanr(input_df, axis=0)[0]
    input_data_uppertria = input_data[np.triu_indices_from(input_data, 1)]
    corrs_df = pd.DataFrame(data=input_data_uppertria,
                            columns=[f'UT'],
                            index=gene_pairs)
    corrs_df.to_csv(prefix_results / celltype / f'{celltype}_correlation.csv')
    # Filter for 0.3 correlation cutoff
    corrs_df = corrs_df[corrs_df['UT'] > 0.3]
    corrs_df.to_csv(prefix_results / celltype / f'{celltype}_correlation_03filtered.csv')
    print(f"Finished {celltype} with time {time() - starttime}")
