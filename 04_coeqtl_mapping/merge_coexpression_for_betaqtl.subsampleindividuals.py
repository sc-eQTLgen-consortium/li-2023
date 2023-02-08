import pandas as pd
from pathlib import Path
import numpy as np
import argparse


def read_numpy(prefix):
    data = np.load(f'{prefix}.npy')
    columns = [item.strip() for item in open(f'{prefix}.cols.txt', 'r').readlines()]
    rows = [item.strip() for item in open(f'{prefix}.rows.txt', 'r').readlines()]
    return pd.DataFrame(data=data, columns=columns, index=rows)


def concat_numpy_files(celltype, condition, res_prefix, num):
    allres = pd.DataFrame()
    for dataset in ['onemillionv2', 'onemillionv3', 'stemiv2', 'ng']:
        if condition =='stimulated' and dataset == 'ng':
            continue
        else:
            numpyfile_path = res_prefix/condition/dataset/f'{condition}_{celltype}.zscores'
            df = read_numpy(numpyfile_path)
            allres = pd.concat([df, allres], axis=1, join='outer')
            print(f'Adding {dataset}, it has shape:', allres.shape)
    allres.sample(num, axis=1).to_csv(res_prefix/condition/f'{condition}_{celltype}.onemillionv23stemiv2ng.{num}randompeople.zscores.tsv.gz',
                  sep='\t', compression='gzip')
    # allres.sample(50).to_csv(res_prefix / condition / f'{condition}_{celltype}.onemillionv23stemiv2ng.50randompeople.zscores.tsv.gz',
    #                           sep='\t', compression='gzip')
    return allres


def argumentsparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--celltype', type=str, dest='celltype')
    parser.add_argument('--condition', type=str, dest='condition')
    parser.add_argument('--num', type=str, dest='num')
    return parser


workdir = Path("./coeqtl_mapping")
res_prefix = workdir/'input/individual_networks/'

args = argumentsparser().parse_args()
celltype, condition, number = args.celltype, args.condition, int(args.num)
_ = concat_numpy_files(celltype, condition, res_prefix, number)