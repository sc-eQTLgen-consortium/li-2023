from pathlib import Path
import numpy as np
import scanpy as sc
import re
import pandas as pd


prefix = Path('./seurat_objects')
data_path_dic = {'onemillionv2':prefix/'1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.sct.h5ad',
                 'stemiv2': prefix / 'cardio.integrated.20210301.stemiv2.h5ad',
                 'onemillionv3': prefix / "1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.SCT.h5ad",
                 'ng': prefix / 'pilot3_seurat3_200420_sct_azimuth.h5ad'}


# extract timepoint from timepoint - stimulation annotation
def get_time(x):
    if x == 'UT':
        return x
    else:
        pattern = re.compile(r'\d+h')
        return re.findall(pattern, x)[0]


def count_nonzeroratio(data_sc):
    df = pd.DataFrame(data=data_sc.X.toarray(),
                      index=data_sc.obs.index,
                      columns=data_sc.var.index)
    nonzerocounts = np.count_nonzero(df.values, axis=0)/df.shape[0]
    return nonzerocounts


def load_onemillion(data_name, data_sc):
    var_df = pd.DataFrame(index=data_sc.var.index.values)
    data_sc.obs['time'] = [get_time(x) for x in data_sc.obs['timepoint']]
    for condition in data_sc.obs['time'].unique():
        for celltype in data_sc.obs['cell_type_lowerres'].unique():
            print(condition, celltype)
            subset_sc = data_sc[(data_sc.obs['time']==condition) &
                                (data_sc.obs['cell_type_lowerres']==celltype)]
            var_df[f'{data_name}_{condition}_{celltype}'] = count_nonzeroratio(subset_sc)
    return var_df


def load_ng(data_sc):
    var_df = pd.DataFrame(index=data_sc.var.index.values)
    celltype_maping = {'CD4 T': 'CD4T', 'CD8 T': 'CD8T', 'Mono': 'monocyte', 'DC': 'DC', 'NK': 'NK',
                       'other T': 'otherT', 'other': 'other', 'B': 'B'}
    data_sc.obs['cell_type_mapped_to_onemillion'] = [celltype_maping.get(name) for name in
                                                     data_sc.obs['predicted.celltype.l1']]
    for celltype in data_sc.obs['cell_type_mapped_to_onemillion'].unique():
        print(celltype)
        subset_sc = data_sc[(data_sc.obs['cell_type_mapped_to_onemillion']==celltype)]
        var_df[f'ng_{celltype}'] = count_nonzeroratio(subset_sc)
    return var_df


def load_stemi(dataname, data_sc):
    var_df = pd.DataFrame(index=data_sc.var.index.values)
    for condition in data_sc.obs['timepoint.final'].unique():
        for celltype in data_sc.obs['cell_type_lowerres'].unique():
            print(condition, celltype)
            subset_sc = data_sc[(data_sc.obs['timepoint.final']==condition) &
                           (data_sc.obs['cell_type_lowerres']==celltype)]
            var_df[f'{dataname}_{condition}_{celltype}'] = count_nonzeroratio(subset_sc)
    return var_df


def get_expressed_ratio(datasetname):
    data_sc = sc.read_h5ad(data_path_dic[datasetname])
    if datasetname.startswith('onemillion'):
        var_df = load_onemillion(datasetname, data_sc)
    elif datasetname.startswith('stemi'):
        var_df = load_stemi(datasetname, data_sc)
    else:
        var_df = load_ng(data_sc)
    return var_df


def calculate_genes_withnonzeroratio(datasetname, savepath):
    print("Processing ", datasetname)
    var_df = get_expressed_ratio(datasetname)
    var_df.to_csv(savepath, sep='\t')
    return var_df


work_dir = Path('/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/')
nonzero_savepath = work_dir/'coeqtl_mapping/input/gene_pair_selection/annotations/'
for datasetname in ['stemiv2', 'ng', 'onemillionv2', 'onemillionv3']:
    print('Processing ', datasetname)
    savepath = nonzero_savepath/f'{datasetname}.genes_nonzeroratio.tsv'
    var_df = calculate_genes_withnonzeroratio(datasetname, savepath)




