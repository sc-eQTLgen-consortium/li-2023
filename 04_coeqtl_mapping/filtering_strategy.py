import pandas as pd
import numpy as np
from pathlib import Path
import argparse
from scipy.stats import norm
from time import time


sig_thres_zscore = norm.ppf(1-0.025)
individual_network_prefix = Path("./input/individual_networks")
saveprefix = Path("./input/gene_pair_selection/annotations/")
def read_numpy(prefix):
    data = np.load(f'{prefix}.npy')
    rows = [item.strip() for item in open(f'{prefix}.rows.txt', 'r').readlines()]
    cols = [item.strip() for item in open(f'{prefix}.cols.txt', 'r').readlines()]
    return pd.DataFrame(data=data, columns=cols, index=rows)


def merge_datasets(celltype, condition):
    res_df = pd.DataFrame()
    for datasetname in ['stemiv2', 'onemillionv2', 'onemillionv3', 'ng']:
        data_path = individual_network_prefix / condition / datasetname / f'{condition}_{celltype}.zscores'
        startime = time()
        df = read_numpy(data_path)
        res_df = pd.concat([res_df, df], axis=1)
        print(f'Merged {datasetname}, it took', time() - startime)
    return res_df


def calculate_significance_freq(zscore_df, thres=sig_thres_zscore):
    freqs = (abs(zscore_df.values) > thres).sum(axis=1)
    assert len(freqs) == zscore_df.shape[0]
    return freqs


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--celltype', dest='celltype')
    parser.add_argument('--condition', dest='condition')
    return parser


def main():
    args = parse().parse_args()
    celltype, condition = args.celltype, args.condition
    celltype_condition_df = merge_datasets(celltype, condition)
    celltype_condition_df['sig_count'] = calculate_significance_freq(celltype_condition_df)
    celltype_condition_df['sig_freq'] = [item/celltype_condition_df.shape[1] for item in celltype_condition_df['sig_count']]
    print(celltype, celltype_condition_df[celltype_condition_df['sig_freq']>=0.1].shape)
    celltype_condition_df[['sig_count', 'sig_freq']].to_csv(saveprefix/f'{condition}_{celltype}.significance_frequency.tsv',
                                               sep='\t')
    return celltype_condition_df


if __name__ == '__main__':
    _ = main()