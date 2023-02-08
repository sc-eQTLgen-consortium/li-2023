####################################################################################
# Calculate per sample correlation with subsampled number of cells per donor
# to explore the relationship between number of cells and
# and concordance between donors
# Individuals with a total number of cells below the respective subsampled value
# are not tested, sampling range from 25 cells to the 75% quantile for the cell type
# (so that at least 25% of the individuals can be included each time)
# Selecting again genes expressed in at least 50% of the cells
# Input: seurat objects with Oelen v2 and v3 dataset
# Output: Csv file with Pearson correlation values for all individual comparison
#         per celltype and subsampled number of cells
####################################################################################

import scanpy as sc
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import pearsonr, spearmanr

def select_gene_nonzeroratio(df, ratio):
    '''
        Select genes with non-zero ratio across all cells > specified ratio
    '''
    nonzerocounts = np.count_nonzero(df.values, axis=0)/df.shape[0]
    selected_genes = df.columns[nonzerocounts>ratio]
    return selected_genes

def calculate_individual_network(individual_df, n_cells=0, random_state=8):
    '''
        Randomly select the n_cells from individual_df to calculate the gene-gene spearman network;
        if n_cell not set, then use all cells from the individual
        Return: A list of correlation coefficients
        '''
    if n_cells > 0:
        specific_individual_network = individual_df.sample(n_cells, random_state=random_state).corr(method='spearman')
    else:
        specific_individual_network = individual_df.corr(method='spearman')
    indices = np.tril_indices_from(np.zeros((specific_individual_network.shape[0],
                                             specific_individual_network.shape[0])),
                                   k=1)
    return specific_individual_network.values[indices].flatten()

def calculate_correlation(celltype_df, celltype_obs, selected_genes, select_individuals, n_cells=100):
    '''
        Calculate all inidividual networks for certain n_cells, and selected_indidviauls
        Input: celltype_df: gene index in columns, cell index in rows, cell index will be used to match the index of
        celltype_df; it needs at least one column named 'assignment', storing the individual index for
        each cell
        selected_genes: a list of gene index
        select_individuals: a list of individual index
        n_cells: number of cells to select from each individual, it should be smaller than the max number of cells
        in all individuals
        Output: all_individuals_correlation: a dataframe, each column of values is all the gene-gene spearman
        correlation coefficients
        correlation_of_individual_correlations: spearman correlation between all pairs of individuals' networks
        '''
    all_individuals_correlation = pd.DataFrame()
    if select_individuals is not None:
        for assignment in select_individuals:
                allcells_individual = celltype_obs[celltype_obs.assignment==assignment].index.values
                specific_individual_df = celltype_df[selected_genes].loc[allcells_individual]
                all_individuals_correlation[assignment] = calculate_individual_network(specific_individual_df, n_cells)
    correlation_of_individual_correlations = all_individuals_correlation.corr(method='pearson')
    individual_indices = np.triu_indices_from(correlation_of_individual_correlations.values, k=1)
    return all_individuals_correlation, correlation_of_individual_correlations.values[individual_indices]

# define path (run for Oelen v2 and v3 dataset separately)
version2 = False
if version2:
    input_path = 'seurat_objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.sct.h5ad'
    output_path ='co-expression_indivs_subsampled/correlation_individuals_subsampled_1M_v2.csv'
else:
    input_path = 'seurat_objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.SCT.h5ad'
    output_path ='co-expression_indivs_subsampled/correlation_individuals_subsampled_1M_v3.csv'

# load single cell data
alldata = sc.read_h5ad(input_path)

# filter to look only at UT cells
alldata = alldata[alldata.obs.timepoint == "UT"].copy()

# select common individual per celltype
celltypes = ['CD4T','NK','monocyte','CD8T','B','DC']
selected_individuals = {}
selected_individuals_cell_number = {}
for celltype in celltypes:
    celltype_data = alldata[alldata.obs.cell_type_lowerres == celltype]
    selected_individuals_cell_number[f'{celltype}'] = celltype_data.obs.assignment.value_counts().values
    selected_individuals[f'{celltype}'] = celltype_data.obs.assignment.value_counts().index
    #Check distribution of cells per individual
    print(celltype_data.obs.assignment.value_counts().describe())

# calculate for each celltype
all_celltype_res = pd.DataFrame()
for celltype in tqdm(celltypes):
    celltype_data = alldata[alldata.obs.cell_type_lowerres == celltype]
    celltype_df = pd.DataFrame(data=celltype_data.X.toarray(),
                               columns=celltype_data.var.index, # genes
                               index=celltype_data.obs.index) # cells
                               
    # select only genes expressed in at least 50% of the cells
    selected_genes = select_gene_nonzeroratio(celltype_df, ratio=0.5)
    
    # Run each cell type for different number of cells so that at least 25% of individuals have that many cells
    for cell_num in range(25, int(np.quantile(selected_individuals_cell_number[f'{celltype}'],0.75)),25):
        
        print(cell_num)
        
        # Select all individuals that have enough cells
        indivs = selected_individuals[f'{celltype}'][selected_individuals_cell_number[f'{celltype}']>=cell_num]
        # Get all pairwise correlation for these pairs
        celltype_correlations = pd.DataFrame(data=calculate_correlation(celltype_df, celltype_data.obs,
                                                                       selected_genes=selected_genes,
                                                                       n_cells=cell_num,
                                                                       select_individuals=indivs)[1],
                                            columns=['corr'])
        celltype_correlations['celltype'] = celltype
        celltype_correlations['cell_num'] = cell_num
        all_celltype_res = pd.concat([all_celltype_res, celltype_correlations], axis=0)

all_celltype_res.to_csv(output_path)
