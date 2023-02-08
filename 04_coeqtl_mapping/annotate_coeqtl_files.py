import pandas as pd
from pathlib import Path
import numpy as np
import argparse


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--celltype', dest='celltype')
    parser.add_argument('--networkcelltype', dest='networkcelltype')
    parser.add_argument('--filtertype', dest='filtertype')
    return parser

args = parse().parse_args()
celltype = args.celltype
filtertype = args.filtertype
networkcelltype = args.networkcelltype
# filtertype = 'filtered_results'
workdir = Path("./coeqtl_mapping/")
coeqtl_filepath = workdir/f'output/{filtertype}/UT_{celltype}/coeqtls_fullresults_fixed.all.tsv.gz'

def find_gene2(genepair, eqtlgene):
    gene1, gene2 = genepair.split(';')
    if gene1 == eqtlgene:
        return gene2
    else:
        return gene1

coeqtl_df = pd.read_csv(coeqtl_filepath, sep='\t', compression='gzip')
coeqtl_df['gene2'] = [find_gene2(item[0], item[1]) for item in coeqtl_df[['Gene', 'eqtlgene']].values]
unique_genepairs = list(set(coeqtl_df['Gene']))


network_prefix = Path("./coeqtl_mapping/input/individual_networks/UT/")
def annotate_by_datasets(datasetname, coeqtl_df, unique_genepairs):
    def read_numpy(prefix):
        data = np.load(f'{prefix}.npy')
        columns = [item.strip() for item in open(f'{prefix}.cols.txt', 'r').readlines()]
        rows = [item.strip() for item in open(f'{prefix}.rows.txt', 'r').readlines()]
        return pd.DataFrame(data=data, columns=columns, index=rows)
    print(f"Loading {datasetname}.")
    network_df = read_numpy(network_prefix / datasetname / f'UT_{networkcelltype}.zscores')
    individual_ids = network_df.columns.copy()
    common_genepairs = list(set(unique_genepairs) & set(network_df.index))
    selected_network_df = network_df.loc[common_genepairs]
    selected_network_df[f'var_{datasetname}'] = np.nanvar(selected_network_df[individual_ids].values, axis=1)
    selected_network_df[f'mean_{datasetname}'] = np.nanmean(selected_network_df[individual_ids].values, axis=1)
    var_mean_dic = selected_network_df[[f'var_{datasetname}', f'mean_{datasetname}']].T.to_dict()
    get_var = lambda x:var_mean_dic.get(x)[f'var_{datasetname}'] if x in var_mean_dic else np.nan
    get_mean = lambda x:var_mean_dic.get(x)[f'mean_{datasetname}'] if x in var_mean_dic else np.nan
    coeqtl_df[f'var_{datasetname}'] = [get_var(genepair) for genepair in coeqtl_df['Gene']]
    coeqtl_df[f'mean_{datasetname}'] = [get_mean(genepair) for genepair in coeqtl_df['Gene']]
    return coeqtl_df

for datasetname in ['onemillionv2', 'onemillionv3', 'stemiv2', 'ng']:
    coeqtl_df = annotate_by_datasets(datasetname, coeqtl_df, unique_genepairs)


def annotate_with_nonzero(df, celltype, datasetname, condition='UT'):
    nonzeroratio_prefix = Path(
        "./coeqtl_mapping/input/gene_pair_selection/annotations/")
    nonzeroratio_path = nonzeroratio_prefix/f'{datasetname}.genes_nonzeroratio.tsv'
    nonzero_df = pd.read_csv(nonzeroratio_path, sep='\t', index_col=0)
    if condition == 'UT' and datasetname == 'stemiv2':
        colname = f'{datasetname}_t8w_{celltype}'
    elif condition == 'UT' and datasetname.startswith('onemillion'):
        colname = f'{datasetname}_UT_{celltype}'
    elif condition == 'UT' and datasetname.startswith('ng'):
        colname = f'{datasetname}_{celltype}'
    else:
        raise NotImplementedError(f"{datasetname} {celltype} not understood")
    nonzero_dict = nonzero_df[colname].T.to_dict()
    df[f'eqtlgene_nonzeroratio_{datasetname}'] = [nonzero_dict.get(genename) for genename in df['eqtlgene']]
    df[f'gene2_nonzeroratio_{datasetname}'] = [nonzero_dict.get(genename) for genename in df['gene2']]
    return df

for datasetname in ['onemillionv2', 'onemillionv3', 'stemiv2', 'ng']:
    print(datasetname)
    coeqtl_df = annotate_with_nonzero(coeqtl_df, networkcelltype, datasetname)


coeqtl_df.to_csv(workdir/f'output/{filtertype}/UT_{celltype}/coeqtls_fullresults_fixed.all.annotated.tsv.gz',
                 compression='gzip', sep='\t', index=False)

