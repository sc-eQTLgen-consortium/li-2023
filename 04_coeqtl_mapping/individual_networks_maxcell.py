import os
import re
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr
from scipy.stats import t, norm
from tqdm import tqdm
import argparse
from scipy.stats import rankdata
from collections import namedtuple


def get_time(x):
    if x == 'UT':
        return x
    else:
        pattern = re.compile(r'\d+h')
        return re.findall(pattern, x)[0]


class DATASET:
    def __init__(self, datasetname):
        self.name = datasetname
        self.path_prefix = Path("./seurat_objects")
        self.information = self.get_information()
    def get_information(self):
        if self.name == 'onemillionv2':
            self.path = '1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.sct.h5ad'
            self.individual_id_col = 'assignment'
            self.timepoint_id_col = 'time'
            self.celltype_id = 'cell_type_lowerres'
            self.chosen_condition = {'UT': 'UT',
                                     'stimulated': '3h'}
        elif self.name == 'onemillionv3':
            self.path = '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.SCT.h5ad'
            self.individual_id_col = 'assignment'
            self.timepoint_id_col = 'time'
            self.celltype_id = 'cell_type_lowerres'
            self.chosen_condition = {'UT': 'UT',
                                     'stimulated': '3h'}
        elif self.name == 'stemiv2':
            self.path = 'cardio.integrated.20210301.stemiv2.h5ad'
            self.individual_id_col = 'assignment.final'
            self.timepoint_id_col = 'timepoint.final'
            self.celltype_id = 'cell_type_lowerres'
            self.chosen_condition = {'UT': 't8w',
                                     'stimulated': 'Baseline'}
        elif self.name == 'ng':
            self.path = 'pilot3_seurat3_200420_sct_azimuth.h5ad'
            self.individual_id_col = 'snumber'
            self.celltype_id = 'cell_type_mapped_to_onemillion'
        else:
            raise IOError("Dataset name not understood.")
    def load_dataset(self):
        self.get_information()
        print(f'Loading dataset {self.name} from {self.path_prefix} {self.path}')
        self.data_sc = sc.read_h5ad(self.path_prefix / self.path)
        if self.name.startswith('onemillion'):
            self.data_sc.obs['time'] = [get_time(item) for item in self.data_sc.obs['timepoint']]
        elif self.name == 'ng':
            celltype_maping = {'CD4 T': 'CD4T', 'CD8 T': 'CD8T', 'Mono': 'monocyte', 'DC': 'DC', 'NK': 'NK',
                               'other T': 'otherT', 'other': 'other', 'B': 'B'}
            self.data_sc.obs['cell_type_mapped_to_onemillion'] = [celltype_maping.get(name) for name in
                                                          self.data_sc.obs['predicted.celltype.l1']]



def select_gene_nonzeroratio(df, ratio):
    nonzerocounts = np.count_nonzero(df.values, axis=0) / df.shape[0]
    selected_genes = df.columns[nonzerocounts > ratio]
    return selected_genes


def corr_to_z(coef, num):
    t_statistic = coef * np.sqrt((num - 2) / (1 - coef ** 2))
    prob = t.cdf(t_statistic, num - 2)
    z_score = norm.ppf(prob)
    positive_coef_probs = 1 - prob
    positive_coef_probs[coef < 0] = 0
    negative_coef_probs = prob
    negative_coef_probs[coef > 0] = 0
    probs = negative_coef_probs + positive_coef_probs
    return z_score, probs


def z_to_corr(z, num):
    prob = norm.cdf(z)
    t_statistic = t.ppf(prob, num - 2)
    corr = t_statistic / np.sqrt(num - 2 + t_statistic ** 2)
    return corr


def get_om_name(filename):
    pattern = re.compile(r'LLDeep_\d\d\d\d')
    return re.findall(pattern, filename)[0]


def get_stemi_name(filename):
    pattern = re.compile(r'TEST_\d.')
    return re.findall(pattern, filename)[0]


def save_numpy(data_df, prefix):
    np.save(f'{prefix}.npy', data_df.values)
    with open(f'{prefix}.cols.txt', 'w') as f:
        f.write('\n'.join(data_df.columns))
    with open(f'{prefix}.rows.txt', 'w') as f:
        f.write('\n'.join(data_df.index))
    return None

def _contains_nan(a, nan_policy='propagate'):
    policies = ['propagate', 'raise', 'omit']
    if nan_policy not in policies:
        raise ValueError("nan_policy must be one of {%s}" %
                         ', '.join("'%s'" % s for s in policies))
    try:
        with np.errstate(invalid='ignore'):
            contains_nan = np.isnan(np.sum(a))
    except TypeError:
        try:
            contains_nan = np.nan in set(a.ravel())
        except TypeError:
            contains_nan = False
            nan_policy = 'omit'
    if contains_nan and nan_policy == 'raise':
        raise ValueError("The input contains nan values")
    return contains_nan, nan_policy


def _chk_asarray(a, axis):
    if axis is None:
        a = np.ravel(a)
        outaxis = 0
    else:
        a = np.asarray(a)
        outaxis = axis
    if a.ndim == 0:
        a = np.atleast_1d(a)
    return a, outaxis


def spearmanr_withnan(a, axis=0, nan_policy='propagate'):
    SpearmanrResult = namedtuple('SpearmanrResult', ('correlation', 'pvalue'))
    if axis is not None and axis > 1:
        raise ValueError("spearmanr only handles 1-D or 2-D arrays, supplied axis argument {}, "
                         "please use only values 0, 1 or None for axis".format(axis))
    a, axisout = _chk_asarray(a, axis)
    if a.ndim > 2:
        raise ValueError("spearmanr only handles 1-D or 2-D arrays")
    n_vars = a.shape[1 - axisout]
    n_obs = a.shape[axisout]
    if n_obs <= 1:
        # Handle empty arrays or single observations.
        return SpearmanrResult(np.nan, np.nan)
    a_contains_nan, nan_policy = _contains_nan(a, nan_policy)
    variable_has_nan = np.zeros(n_vars, dtype=bool)
    if a_contains_nan:
        if nan_policy == 'propagate':
            if a.ndim == 1 or n_vars <= 2:
                return SpearmanrResult(np.nan, np.nan)
            else:
                variable_has_nan = np.isnan(a).sum(axis=axisout)
    a_ranked = np.apply_along_axis(rankdata, axisout, a)
    rs = np.corrcoef(a_ranked, rowvar=axisout)
    dof = n_obs - 2  # degrees of freedom
    # rs can have elements equal to 1, so avoid zero division warnings
    with np.errstate(divide='ignore'):
        t_ = rs * np.sqrt((dof/((rs+1.0)*(1.0-rs))).clip(0))
    prob = 2 * t.sf(np.abs(t_), dof)
    # For backwards compatibility, return scalars when comparing 2 columns
    if rs.shape == (2, 2):
        return SpearmanrResult(rs[1, 0], prob[1, 0])
    else:
        rs[variable_has_nan, :] = np.nan
        rs[:, variable_has_nan] = np.nan
        return SpearmanrResult(rs, prob)

def read_numpy(prefix):
    data = np.load(f'{prefix}.npy')
    columns = [item.strip() for item in open(f'{prefix}.rows.txt', 'r').readlines()]
    return pd.DataFrame(data=data, columns=columns, index=columns)


def read_all_files(prefix, genepairs):
    res_df = pd.DataFrame(index=genepairs)
    for filename in os.listdir(prefix):
        if filename.endswith('_coefs.npy'):
            data = np.load(f'{prefix}/{filename}')
            if len(data.shape) > 1:
                data_uppertria = data[np.triu_indices_from(data, 1)]
                individual_id = get_stemi_name(filename)
                res_df[individual_id] = data_uppertria
    return res_df


def get_unique_genepairs(genepair_list, sep=';'):
    unique_pairs = set()
    for genepair in genepair_list:
        reverse_genepair = sep.join(genepair.split(sep))
        if genepair in unique_pairs or reverse_genepair in unique_pairs:
            continue
        else:
            unique_pairs.add(genepair)
    return unique_pairs


def get_genes(genepair_list, sep=';'):
    genes = list(set([gene for genepair in genepair_list for gene in genepair.split(sep)]))
    return genes


def get_genepairs(genelist_path):
    genelist = [item.strip() for item in open(genelist_path, 'r').readlines()]
    genepairs = [';'.join(sorted(item)) for item in combinations(genelist, 2)]
    return genelist, genepairs


def get_individual_networks_selected_genepairs(data_sc, individual_colname, selected_genepairs, maxcell):
    data_df = pd.DataFrame(data=data_sc.X.toarray(),
                           index=data_sc.obs.index,
                           columns=data_sc.var.index)
    selected_genes = list(set([ele for item in selected_genepairs
                               for ele in item.split(';')]) & set(data_sc.var.index))
    selected_genes_sorted_genepairs = [';'.join(sorted(item)) for item in combinations(selected_genes, 2)]
    common_genepairs = list(set(selected_genes_sorted_genepairs) & set(selected_genepairs))
    coef_df = pd.DataFrame(index=common_genepairs)
    coef_p_df = pd.DataFrame(index=common_genepairs)
    zscore_df = pd.DataFrame(index=common_genepairs)
    zscore_p_df = pd.DataFrame(index=common_genepairs)
    data_selected_df = data_df[selected_genes]
    print(f"Begin calculating networks for {len(data_sc.obs[individual_colname].unique())} individuals.")
    for ind_id in tqdm(data_sc.obs[individual_colname].unique()):
        cell_num = data_sc.obs[data_sc.obs[individual_colname] == ind_id].shape[0]
        if cell_num > 10:
            if maxcell>0 and cell_num >= maxcell:
                individual_df = data_selected_df.loc[data_sc.obs[individual_colname] == ind_id].sample(maxcell, random_state=5)
                cell_num = maxcell
            else:
                individual_df = data_selected_df.loc[data_sc.obs[individual_colname] == ind_id]
            # individual_df = data_selected_df.loc[data_sc.obs[individual_colname] == ind_id]
            individual_coefs, individual_coef_ps = spearmanr_withnan(individual_df.values, axis=0)
            try:
                individual_coefs_flatten = pd.DataFrame(data=individual_coefs[np.triu_indices_from(individual_coefs, 1)],
                                                        index=selected_genes_sorted_genepairs).loc[common_genepairs]
                individual_coef_ps_flatten = pd.DataFrame(data=individual_coef_ps[np.triu_indices_from(individual_coefs, 1)],
                                                          index=selected_genes_sorted_genepairs).loc[common_genepairs]
                individual_zscores_flatten, individual_zscore_ps_flatten = corr_to_z(individual_coefs_flatten.values,
                                                                                     cell_num)
                coef_df[ind_id] = individual_coefs_flatten
                coef_p_df[ind_id] = individual_coef_ps_flatten
                zscore_df[ind_id] = individual_zscores_flatten
                zscore_p_df[ind_id] = individual_zscore_ps_flatten
            except:
                continue
        else:
            print("Deleted this individual because of low cell number", cell_num)
    return coef_df, coef_p_df, zscore_df, zscore_p_df


def get_individual_networks_given_celltype_condition_datasetname(celltype, datasetname, condition='UT', maxcell=-1):
    # load the data and data information
    dataset = DATASET(datasetname)
    dataset.load_dataset()
    print(f"{datasetname} loaded.")
    # calculate the individual network for specific condition and celltype
    print(datasetname, celltype, condition)
    work_prefix = Path('./')
    selected_genepairs_path = work_prefix / f'coeqtl_mapping/input/snp_genepair_selection/{condition}_{celltype}_{datasetname}.baseline.tsv'
    selected_genepairs = pd.read_csv(selected_genepairs_path, sep='\t')['genepair_sorted'].values
    if datasetname == 'ng':
        data_selected = dataset.data_sc[(dataset.data_sc.obs[dataset.celltype_id] == celltype)]
    else:
        data_selected = dataset.data_sc[(dataset.data_sc.obs[dataset.celltype_id] == celltype) &
                                        (dataset.data_sc.obs[dataset.timepoint_id_col] == dataset.chosen_condition[condition])]
    individual_coefs_df, individual_coefs_p_df, individual_zscores_df, individual_zscores_p_df = \
        get_individual_networks_selected_genepairs(
        data_selected,
        dataset.individual_id_col,
        selected_genepairs,
        maxcell
    )
    print(individual_coefs_df.head())
    save_prefix = Path('./coeqtl_mapping/input')
    if not os.path.exists(save_prefix / 'individual_networks' / condition / datasetname):
        os.mkdir(save_prefix / 'individual_networks' / condition / datasetname)
    save_numpy(individual_zscores_df,
               save_prefix / 'individual_networks' /  condition / datasetname / f'{condition}_{celltype}.max{maxcell}cells.zscores')
    print("Saved ")
    return individual_coefs_df, individual_coefs_p_df, individual_zscores_df, individual_zscores_p_df


def argumentsparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--datasetname', type=str, dest='datasetname')
    parser.add_argument('--celltype', type=str, dest='celltype')
    parser.add_argument('--condition', type=str, dest='condition')
    parser.add_argument('--maxcell', type=float, dest='maxcell')
    return parser

def run_get_individual_networks_given_celltype_condition_datasetname():
    args = argumentsparser().parse_args()
    print(f"Starting to calculate individual network for {args.datasetname}, {args.celltype}, {args.condition}, "
          f"for max cell number {args.maxcell}.")
    _ = get_individual_networks_given_celltype_condition_datasetname(celltype=args.celltype,
                                                                     condition=args.condition,
                                                                     datasetname=args.datasetname,
                                                                     maxcell=int(args.maxcell))
    return None

if __name__ == '__main__':
    run_get_individual_networks_given_celltype_condition_datasetname()