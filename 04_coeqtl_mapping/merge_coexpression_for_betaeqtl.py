import pandas as pd
from pathlib import Path
import numpy as np
import argparse


def read_numpy(prefix):
    data = np.load(f'{prefix}.npy')
    columns = [item.strip() for item in open(f'{prefix}.cols.txt', 'r').readlines()]
    rows = [item.strip() for item in open(f'{prefix}.rows.txt', 'r').readlines()]
    return pd.DataFrame(data=data, columns=columns, index=rows)


def concat_numpy_files(celltype, condition, res_prefix):
    allres = pd.DataFrame()
    for dataset in ['onemillionv2', 'onemillionv3', 'stemiv2', 'ng']:
        if condition =='stimulated' and dataset == 'ng':
            continue
        else:
            numpyfile_path = res_prefix/condition/dataset/f'{condition}_{celltype}.zscores'
            df = read_numpy(numpyfile_path)
            allres = pd.concat([df, allres], axis=1)
            print(f'Adding {dataset}, it has shape:', allres.shape)
    allres.to_csv(res_prefix/condition/f'{condition}_{celltype}.onemillionv23stemiv2ng.zscores.tsv', sep='\t')
    return allres


def argumentsparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--celltype', type=str, dest='celltype')
    parser.add_argument('--condition', type=str, dest='condition')
    return parser


workdir = Path("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping")
res_prefix = workdir/'input/individual_networks/'

args = argumentsparser().parse_args()
celltype, condition = args.celltype, args.condition
_ = concat_numpy_files(celltype, condition, res_prefix)