import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm
import gzip


workdir = Path('/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping')
annotation_path = '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/probeannotation/singleCell-annotation-stripped.tsv'
mappingdic = pd.read_csv('/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/resources/features_v3_reformated_names.tsv',
                          sep='\t', names=['geneid', 'genename']).set_index('geneid')['genename'].T.to_dict()
annotation_df = pd.read_csv(annotation_path, sep='\t')
annotation_df['chr_pos'] = ['_'.join([str(ele) for ele in item]) for item in annotation_df[['Chr', 'ChrStart', 'ChrEnd']].values]
annotation_df['genename'] = [mappingdic.get(ensemblid) for ensemblid in annotation_df['Ensembl']]
annotation_dict = annotation_df.set_index('chr_pos')['genename'].T.to_dict()

def update_perm(old_p_list, new_p_list):
    return np.min([old_p_list, new_p_list], axis=0)


def find_eqtlsnp_gene(snp, genepair, coeqtl_annotation_dic):
    genepair_chrpos = coeqtl_annotation_dic.get(genepair)
    eqtlgene = annotation_dict.get(genepair_chrpos)
    snp_genepair = '_'.join([snp, eqtlgene])
    return snp_genepair


def loop_through_one_batch_perm(batch_perm_path, snpgene1_minpvalues_dict, coeqtl_annotation_dict):
    with gzip.open(batch_perm_path, 'rb') as f:
        f.readline()
        while True:
            line = f.readline().decode('utf-8')
            if not line:
                break
            else:
                linecontent = line.strip().split('\t')
                perm_ps = [float(ele) for ele in linecontent[2:102]]
                snp_gene1 = find_eqtlsnp_gene(linecontent[1], linecontent[0], coeqtl_annotation_dict)
                snpgene1_minpvalues_dict[snp_gene1] = update_perm(snpgene1_minpvalues_dict[snp_gene1],
                                                                      perm_ps)
    return snpgene1_minpvalues_dict


def update_dictionary_per_permutation_batch(batch_perm_path, snpgene1_minpvalues_df, coeqtl_annotation_dict):
    batch_perm_df = pd.read_csv(batch_perm_path, compression='gzip', sep='\t')
    # print(batch_perm_df.head())
    batch_perm_df['chr_pos'] = [coeqtl_annotation_dict.get(genepair) for genepair in batch_perm_df['Gene']]
    batch_perm_df['eqtlgene'] = [annotation_dict.get(chrpos) for chrpos in batch_perm_df['chr_pos']]
    batch_perm_df['snp_eqtlgene'] = ['_'.join(item) for item in batch_perm_df[['SNP', 'eqtlgene']].values]
    # print(batch_perm_df.head())
    merge_columns = ['snp_eqtlgene'] + [f'Perm{ind}' for ind in range(100)]
    merged_df = pd.concat([batch_perm_df[merge_columns], snpgene1_minpvalues_df[merge_columns]],
                          axis=0)
    # print(merged_df.head())
    reduced_df = merged_df.groupby(by='snp_eqtlgene').agg(min)
    reduced_df['snp_eqtlgene'] = reduced_df.index
    # print(reduced_df.head())
    return reduced_df


def save_numpy(data_df, prefix):
    np.save(f'{prefix}.npy', data_df.values)
    with open(f'{prefix}.cols.txt', 'w') as f:
        f.write('\n'.join([str(ele) for ele in data_df.columns]))
    with open(f'{prefix}.rows.txt', 'w') as f:
        f.write('\n'.join([str(ele) for ele in data_df.index]))
    return None


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--eqtl_path', dest='eqtl_path')
    parser.add_argument('--result_prefix', dest='result_prefix')
    parser.add_argument('--save_prefix', dest='save_prefix')
    parser.add_argument('--annotation_prefix', dest='annotation_prefix')
    return parser



def main():
    args = arguments().parse_args()
    eqtl_path, results_prefix, save_prefix = args.eqtl_path, Path(args.result_prefix), Path(args.save_prefix)
    # load eqtl path
    eqtl_df = pd.read_csv(eqtl_path, sep='\t')
    eqtl_df['snp_gene1'] = ['_'.join(item) for item in eqtl_df[['SNPName', 'genename']].values]
    unique_snpgene1 = eqtl_df['snp_gene1'].values
    # initialize the dict to contain the
    # snpgene1_minpvalues_dict = {item: np.ones(100) for item in unique_snpgene1}
    snpgene1_minpvalues_df = pd.DataFrame(data=np.ones((len(unique_snpgene1), 100)),
                                          columns=[f'Perm{ind}' for ind in range(100)])
    snpgene1_minpvalues_df['snp_eqtlgene'] = unique_snpgene1
    # loop through all batch permutation files
    coeqtl_annotation_path = f'{args.annotation_prefix}.genepairs.annotation.gene1position.noduplicated.tsv'
    coeqtl_annotation_df = pd.read_csv(coeqtl_annotation_path, sep='\t')
    coeqtl_annotation_df['chr_pos'] = ['_'.join([str(ele) for ele in item]) for item in
                                       coeqtl_annotation_df[['Chr', 'ChrStart', 'ChrEnd']].values]
    coeqtl_annotation_dict = coeqtl_annotation_df.set_index('ArrayAddress')['chr_pos'].T.to_dict()
    for filename in tqdm(os.listdir(results_prefix / 'noduplicated/output')):
        if '-Permutations.txt.gz' in filename:
            snpgene1_minpvalues_df = update_dictionary_per_permutation_batch(results_prefix / 'noduplicated/output'/filename,
                                                                       snpgene1_minpvalues_df, coeqtl_annotation_dict)
    coeqtl_annotation_path = f'{args.annotation_prefix}.genepairs.annotation.gene1position.duplicatedversion1.tsv'
    coeqtl_annotation_df = pd.read_csv(coeqtl_annotation_path, sep='\t')
    coeqtl_annotation_df['chr_pos'] = ['_'.join([str(ele) for ele in item]) for item in
                                       coeqtl_annotation_df[['Chr', 'ChrStart', 'ChrEnd']].values]
    coeqtl_annotation_dict = coeqtl_annotation_df.set_index('ArrayAddress')['chr_pos'].T.to_dict()
    for filename in tqdm(os.listdir(results_prefix / 'duplicatedversion1/output')):
        if '-Permutations.txt.gz' in filename:
            snpgene1_minpvalues_df = update_dictionary_per_permutation_batch(results_prefix / 'duplicatedversion1/output'/filename,
                                                                   snpgene1_minpvalues_df, coeqtl_annotation_dict)
    coeqtl_annotation_path = f'{args.annotation_prefix}.genepairs.annotation.gene1position.duplicatedversion2.tsv'
    coeqtl_annotation_df = pd.read_csv(coeqtl_annotation_path, sep='\t')
    coeqtl_annotation_df['chr_pos'] = ['_'.join([str(ele) for ele in item]) for item in
                                       coeqtl_annotation_df[['Chr', 'ChrStart', 'ChrEnd']].values]
    coeqtl_annotation_dict = coeqtl_annotation_df.set_index('ArrayAddress')['chr_pos'].T.to_dict()
    for filename in tqdm(os.listdir(results_prefix / 'duplicatedversion2/output')):
        if '-Permutations.txt.gz' in filename:
            snpgene1_minpvalues_df = update_dictionary_per_permutation_batch(results_prefix / 'duplicatedversion2/output'/filename,
                                                                   snpgene1_minpvalues_df, coeqtl_annotation_dict)
    # snpgene1_minpvalues_df = pd.DataFrame.from_dict(snpgene1_minpvalues_dict)
    snpgene1_minpvalues_df.to_csv(save_prefix / 'concated_alltests_permutations_fixed.tsv.gz',
                                    sep='\t', compression='gzip')
    return snpgene1_minpvalues_df


if __name__ == '__main__':
    _ = main()
