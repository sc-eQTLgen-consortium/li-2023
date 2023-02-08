import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr
import argparse


def argumentsparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--filtertype', type=str, dest='filtertype')
    return parser

def prepare_for_rb_BIOS_replication(celltype, filtertype, bios_replication_type='onlyRNAAlignMetrics'):
    '''
    Rb Calculation preparation for BIOS replication
    '''
    workdir = Path("./coeqtl_mapping")
    coeqtl_path = workdir/f'output/{filtertype}/UT_{celltype}/coeqtls_fullresults_fixed.sig.withbios{bios_replication_type}.tsv.gz'
    coeqtl_df = pd.read_csv(coeqtl_path, sep='\t', compression='gzip')
    coeqtl_df['theta'] = 0
    def flip_direction(allele1, allele2, coef2):
        if allele1 == allele2:
            return coef2
        else:
            return -1*coef2
    coeqtl_df['flipped_bios_beta'] = [flip_direction(item[0],
                                                     item[1],
                                                     item[2]) for item in
                                      coeqtl_df[['SNPEffectAllele',
                                                 'assessed_allele_bios',
                                                 'coef_bios']].values]
    coeqtl_df[['snp_genepair', 'snp_eqtlgene',
               'flipped_bios_beta', 'std err_bios',
               'MetaBeta', 'MetaSE', 'theta']].dropna().to_csv(
        workdir/f'bios/{bios_replication_type}/{filtertype}/UT_{celltype}/replication_parameters.csv'
    )
    return coeqtl_df


def find_gene2(genepair, eqtlgene):
    gene1, gene2 = genepair.split(';')
    if gene1 == eqtlgene:
        return gene2
    else:
        return gene1


def flip_direction(df, flipcol, allele1_col, allele2_col):
    df = df.rename({flipcol: f'{flipcol}_ori'}, axis=1)
    def flip(x1, x2, x3):
        if not pd.isnull(x1):
            if x2 == x3:
                return x1
            else:
                return -x1
        else:
            return x1
    df[f'{flipcol}'] = [flip(score, allele1, allele2) for (score, allele1,allele2)
                        in df[[f'{flipcol}_ori', allele1_col, allele2_col]].values]
    return df


# coeQTLs
args = argumentsparser().parse_args()
filtertype = args.filtertype
workdir = Path("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping")
celltypes = ['CD4T', 'CD8T', 'monocyte', 'DC', 'B', 'NK']
for celltype_replication in celltypes:
    print(f"Discovery: {celltype_replication}")
    replication_coeqtl_path = workdir / f'output/{filtertype}/UT_{celltype_replication}/coeqtls_fullresults_fixed.all.tsv.gz'
    replication_coeqtl_df = pd.read_csv(replication_coeqtl_path, sep='\t', compression='gzip')
    replication_coeqtl_df['gene2'] = [find_gene2(x[0], x[1]) for x in
                                      replication_coeqtl_df[['Gene',
                                                             'eqtlgene']].values]
    replication_coeqtl_df['snp_eqtlgene_gene2'] = ['_'.join([item[0], item[1]]) for item in
                                                replication_coeqtl_df[['snp_eqtlgene',
                                                                       'gene2']].values]
    replication_coeqtl_df = replication_coeqtl_df.set_index('snp_eqtlgene_gene2')
    replication_coexpression_df = pd.read_csv(workdir/f'input/individual_networks/UT/UT_{celltype_replication}.sigcoeQTLs.tsv.gz',
                         compression='gzip', sep='\t', index_col=0)
    for celltype_discovery in celltypes:
        if celltype_replication != celltype_discovery:
            print(f"Replication: {celltype_discovery}")
            discovery_coeqtl_path = workdir / f'output/{filtertype}/UT_{celltype_discovery}/coeqtls_fullresults_fixed.sig.tsv.gz'
            discovery_coeqtl_df = pd.read_csv(discovery_coeqtl_path, sep='\t', compression='gzip')
            discovery_coeqtl_df['gene2'] = [find_gene2(x[0], x[1]) for x in
                                              discovery_coeqtl_df[['Gene',
                                                                     'eqtlgene']].values]
            discovery_coeqtl_df['snp_eqtlgene_gene2'] = ['_'.join([item[0], item[1]]) for item in
                                                           discovery_coeqtl_df[['snp_eqtlgene',
                                                                                  'gene2']].values]
            discovery_coeqtl_df = discovery_coeqtl_df.set_index('snp_eqtlgene_gene2')
            tested_coeqtls = list(set(replication_coeqtl_df.index) & set(discovery_coeqtl_df.index))
            merged_coeqtl_df = pd.concat([replication_coeqtl_df.loc[tested_coeqtls],
                                          discovery_coeqtl_df.loc[tested_coeqtls].add_suffix('_replication')], # todo: here is wrong.. should be discovery
                                         axis=1)
            merged_coeqtl_df = flip_direction(merged_coeqtl_df,
                                              'MetaBeta_replication',
                                              'SNPEffectAllele',
                                              'SNPEffectAllele_replication') # MetaBeta, MetaSE, MetaBeta_replication, MetaSE_replication
            disovery_coexpression_df = pd.read_csv(
                workdir / f'input/individual_networks/UT/UT_{celltype_discovery}.sigcoeQTLs.tsv.gz',
                compression='gzip', sep='\t', index_col=0)
            # find overlapping individuals
            tested_genepairs = list(merged_coeqtl_df['Gene'].unique())
            tested_coexpression_discovery_df = disovery_coexpression_df.loc[tested_genepairs]
            tested_coexpression_discovery_df.replace([np.inf, -np.inf], np.nan, inplace=True)
            tested_coexpression_replication_df = replication_coexpression_df.loc[tested_genepairs]
            tested_coexpression_replication_df.replace([np.inf, -np.inf], np.nan, inplace=True)
            other_col_dict = {genepair:np.nan for genepair in tested_genepairs}
            for genepair in tested_genepairs:
                tested_coexpression_discovery_genepair_nonan = tested_coexpression_discovery_df.loc[genepair].dropna()
                tested_coexpression_replication_genepair_nonan = tested_coexpression_replication_df.loc[genepair].dropna()
                common_individuals = list(set(tested_coexpression_discovery_genepair_nonan.index) & set(tested_coexpression_replication_genepair_nonan.index))
                num_common = len(common_individuals)
                num_discovery = tested_coexpression_discovery_genepair_nonan.shape[0]
                num_replication = tested_coexpression_replication_genepair_nonan.shape[0]
                rho = pearsonr(tested_coexpression_discovery_genepair_nonan[common_individuals],
                               tested_coexpression_replication_genepair_nonan[common_individuals])[0]
                other_col_dict[genepair] = rho * num_common / np.sqrt(num_discovery * num_replication)
            merged_coeqtl_df['theta'] = [other_col_dict.get(genepair) for genepair in merged_coeqtl_df['Gene']]
            merged_coeqtl_df[['MetaBeta',
                              'MetaBeta_replication',
                              'MetaSE',
                              'MetaSE_replication',
                              'theta']].to_csv(workdir/f'output/{filtertype}/rb_calculations/discovery_{celltype_discovery}_replication_{celltype_replication}.tsv.gz',
                                               sep='\t',
                                               compression='gzip')
        else:
            continue


# cmono ncmono and monocyte
filtertype = 'filtered_results'
workdir = Path("./coeqtl_mapping")
celltypes = ['monocyte', 'cMono', 'ncMono']
for celltype_replication in celltypes:
    print(f"Discovery: {celltype_replication}")
    replication_coeqtl_path = workdir / f'output/{filtertype}/UT_{celltype_replication}/coeqtls_fullresults_fixed.all.tsv.gz'
    replication_coeqtl_df = pd.read_csv(replication_coeqtl_path, sep='\t', compression='gzip')
    replication_coeqtl_df['gene2'] = [find_gene2(x[0], x[1]) for x in
                                      replication_coeqtl_df[['Gene',
                                                             'eqtlgene']].values]
    replication_coeqtl_df['snp_eqtlgene_gene2'] = ['_'.join([item[0], item[1]]) for item in
                                                replication_coeqtl_df[['snp_eqtlgene',
                                                                       'gene2']].values]
    replication_coeqtl_df = replication_coeqtl_df.set_index('snp_eqtlgene_gene2')
    replication_coexpression_df = pd.read_csv(workdir/f'input/individual_networks/UT/monocyte_subcelltypes/UT_{celltype_replication}.sigcoeQTLs.tsv.gz',
                         compression='gzip', sep='\t', index_col=0)
    for celltype_discovery in celltypes:
        if celltype_replication != celltype_discovery:
            print(f"Replication: {celltype_discovery}")
            discovery_coeqtl_path = workdir / f'output/{filtertype}/UT_{celltype_discovery}/coeqtls_fullresults_fixed.sig.tsv.gz'
            discovery_coeqtl_df = pd.read_csv(discovery_coeqtl_path, sep='\t', compression='gzip')
            discovery_coeqtl_df['gene2'] = [find_gene2(x[0], x[1]) for x in
                                              discovery_coeqtl_df[['Gene',
                                                                     'eqtlgene']].values]
            discovery_coeqtl_df['snp_eqtlgene_gene2'] = ['_'.join([item[0], item[1]]) for item in
                                                           discovery_coeqtl_df[['snp_eqtlgene',
                                                                                  'gene2']].values]
            discovery_coeqtl_df = discovery_coeqtl_df.set_index('snp_eqtlgene_gene2')
            tested_coeqtls = list(set(replication_coeqtl_df.index) & set(discovery_coeqtl_df.index))
            merged_coeqtl_df = pd.concat([replication_coeqtl_df.loc[tested_coeqtls],
                                          discovery_coeqtl_df.loc[tested_coeqtls].add_suffix('_replication')], # todo: also here it is wrong...
                                         axis=1)
            merged_coeqtl_df = flip_direction(merged_coeqtl_df,
                                              'MetaBeta_replication',
                                              'SNPEffectAllele',
                                              'SNPEffectAllele_replication') # MetaBeta, MetaSE, MetaBeta_replication, MetaSE_replication
            disovery_coexpression_df = pd.read_csv(
                workdir / f'input/individual_networks/UT/monocyte_subcelltypes/UT_{celltype_discovery}.sigcoeQTLs.tsv.gz',
                compression='gzip', sep='\t', index_col=0)
            # find overlapping individuals
            tested_genepairs = list(merged_coeqtl_df['Gene'].unique())
            tested_coexpression_discovery_df = disovery_coexpression_df.loc[tested_genepairs]
            tested_coexpression_discovery_df.replace([np.inf, -np.inf], np.nan, inplace=True)
            tested_coexpression_replication_df = replication_coexpression_df.loc[tested_genepairs]
            tested_coexpression_replication_df.replace([np.inf, -np.inf], np.nan, inplace=True)
            other_col_dict = {genepair:np.nan for genepair in tested_genepairs}
            for genepair in tested_genepairs:
                tested_coexpression_discovery_genepair_nonan = tested_coexpression_discovery_df.loc[genepair].dropna()
                tested_coexpression_replication_genepair_nonan = tested_coexpression_replication_df.loc[genepair].dropna()
                common_individuals = list(set(tested_coexpression_discovery_genepair_nonan.index) & set(tested_coexpression_replication_genepair_nonan.index))
                num_common = len(common_individuals)
                num_discovery = tested_coexpression_discovery_genepair_nonan.shape[0]
                num_replication = tested_coexpression_replication_genepair_nonan.shape[0]
                rho = pearsonr(tested_coexpression_discovery_genepair_nonan[common_individuals],
                               tested_coexpression_replication_genepair_nonan[common_individuals])[0]
                other_col_dict[genepair] = rho * num_common / np.sqrt(num_discovery * num_replication)
            merged_coeqtl_df['theta'] = [other_col_dict.get(genepair) for genepair in merged_coeqtl_df['Gene']]
            merged_coeqtl_df[['MetaBeta',
                              'MetaBeta_replication',
                              'MetaSE',
                              'MetaSE_replication',
                              'theta']].to_csv(workdir/f'output/{filtertype}/rb_calculations/monocyte_subcelltypes/discovery_{celltype_discovery}_replication_{celltype_replication}.tsv.gz',
                                               sep='\t',
                                               compression='gzip')
        else:
            continue


# eQTLs
workdir = Path("./cis_eqtl_single_cell/EMP_mapping_1_12_2021_perm1000/output/")
celltypes = ['CD4T', 'CD8T', 'monocyte', 'DC', 'B', 'NK']
genename_dict = pd.read_csv('/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/resources/features_v3_reformated_names.tsv',
                            sep='\t', names=['ensemblid', 'genename']).set_index('ensemblid')['genename'].T.to_dict()

def read_alldataset_celltype_expression_df(celltype):
    expression_prefix = Path('./expression_files/sources_for_coeqtl')
    df = pd.DataFrame()
    for datasetname in ['1m_v2', '1m_v3', 'NG', 't8w']:
        dataset_df = pd.read_csv(expression_prefix/f'{datasetname}/{celltype}_expression.tsv', sep='\t', index_col=0)
        dataset_df['genename'] = [genename_dict.get(geneid) for geneid in dataset_df.index]
        dataset_df = dataset_df.dropna(subset=['genename']).set_index('genename')
        df = pd.concat([df, dataset_df], axis=1)
    return df

for celltype_replication in celltypes:
    print(f"Discovery: {celltype_replication}")
    replication_coeqtl_path = workdir / f'{celltype_replication}/eQTLsFDR-ProbeLevel.txt.gz'
    replication_coeqtl_df = pd.read_csv(replication_coeqtl_path, sep='\t', compression='gzip')
    replication_coeqtl_df['genename'] = [genename_dict.get(ensemblid) for ensemblid in replication_coeqtl_df['ProbeName']]
    replication_coeqtl_df['snp_gene'] = ['_'.join(item) for item in replication_coeqtl_df[['SNPName', 'genename']].values]
    replication_coeqtl_df = replication_coeqtl_df.set_index('snp_gene')
    replication_coeqtl_df['metabeta'] = [float(item.split(' (')[0]) for item in replication_coeqtl_df['Meta-Beta (SE)']]
    replication_coeqtl_df['SE'] = [float(item.split(' (')[1][:-2]) for item in replication_coeqtl_df['Meta-Beta (SE)']]
    replication_coexpression_df = read_alldataset_celltype_expression_df(celltype_replication)
    for celltype_discovery in celltypes:
        if celltype_replication != celltype_discovery:
            print(f"Replication: {celltype_discovery}")
            discovery_coeqtl_path = workdir / f'{celltype_discovery}/eQTLsFDR0.05-ProbeLevel.txt.gz'
            discovery_coeqtl_df = pd.read_csv(discovery_coeqtl_path, sep='\t', compression='gzip')
            discovery_coeqtl_df['genename'] = [genename_dict.get(ensemblid) for ensemblid in
                                                 discovery_coeqtl_df['ProbeName']]
            discovery_coeqtl_df['snp_gene'] = ['_'.join(item) for item in
                                                 discovery_coeqtl_df[['SNPName', 'genename']].values]
            discovery_coeqtl_df = discovery_coeqtl_df.set_index('snp_gene')
            discovery_coeqtl_df['metabeta'] = [float(item.split(' (')[0]) for item in discovery_coeqtl_df['Meta-Beta (SE)']]
            discovery_coeqtl_df['SE'] = [float(item.split(' (')[1][:-2]) for item in discovery_coeqtl_df['Meta-Beta (SE)']]
            tested_eqtls = list(set(replication_coeqtl_df.index) & set(discovery_coeqtl_df.index))
            merged_coeqtl_df = pd.concat([replication_coeqtl_df.loc[tested_eqtls],
                                          discovery_coeqtl_df.loc[tested_eqtls].add_suffix('_replication')], # todo here it is wrong...
                                         axis=1)
            merged_coeqtl_df = flip_direction(merged_coeqtl_df,
                                              'metabeta_replication',
                                              'AlleleAssessed',
                                              'AlleleAssessed_replication')
            discovery_coexpression_df = read_alldataset_celltype_expression_df(celltype_discovery)
            # find overlapping individuals
            tested_genepairs = list(merged_coeqtl_df['genename'].unique())
            tested_coexpression_discovery_df = discovery_coexpression_df.loc[tested_genepairs]
            tested_coexpression_discovery_df.replace([np.inf, -np.inf], np.nan, inplace=True)
            tested_coexpression_replication_df = replication_coexpression_df.loc[tested_genepairs]
            tested_coexpression_replication_df.replace([np.inf, -np.inf], np.nan, inplace=True)
            other_col_dict = {genepair:np.nan for genepair in tested_genepairs}
            for genepair in tested_genepairs:
                tested_coexpression_discovery_genepair_nonan = tested_coexpression_discovery_df.loc[genepair].dropna()
                tested_coexpression_replication_genepair_nonan = tested_coexpression_replication_df.loc[genepair].dropna()
                common_individuals = list(set(tested_coexpression_discovery_genepair_nonan.index) & set(tested_coexpression_replication_genepair_nonan.index))
                num_common = len(common_individuals)
                num_discovery = tested_coexpression_discovery_genepair_nonan.shape[0]
                num_replication = tested_coexpression_replication_genepair_nonan.shape[0]
                rho = pearsonr(tested_coexpression_discovery_genepair_nonan[common_individuals],
                               tested_coexpression_replication_genepair_nonan[common_individuals])[0]
                other_col_dict[genepair] = rho * num_common / np.sqrt(num_discovery * num_replication)
            merged_coeqtl_df['theta'] = [other_col_dict.get(genepair) for genepair in merged_coeqtl_df['genename']]
            merged_coeqtl_df[['metabeta',
                              'metabeta_replication',
                              'SE',
                              'SE_replication',
                              'theta']].to_csv(f'./coeqtl_mapping/input/snp_selection/rb_calculations/discovery_{celltype_discovery}_replication_{celltype_replication}.tsv.gz',
                                               sep='\t',
                                               compression='gzip')
        else:
            continue



filtertype = 'filtered_results'
workdir = Path("./coeqtl_mapping")
celltypes = ['CD4T', 'CD8T', 'monocyte', 'DC', 'B', 'NK']
for celltype_replication in celltypes:
    for celltype_discovery in celltypes:
        if celltype_replication != celltype_discovery:
            print(celltype_discovery, celltype_replication)
            merged_coeqtl_df = pd.read_csv(workdir/f'output/{filtertype}/rb_calculations/discovery_{celltype_discovery}_replication_{celltype_replication}.tsv.gz',
                                                   sep='\t',
                                                   compression='gzip')
            merged_coeqtl_df = merged_coeqtl_df.rename({
                'MetaBeta_replication': 'MetaBeta_discovery',
                'MetaSE_replication': 'MetaSE_discovery'
            },
            axis=1)
            merged_coeqtl_df = merged_coeqtl_df.rename({
                'MetaBeta': 'MetaBeta_replication',
                'MetaSE': 'MetaSE_replication'
            },
                axis=1)
            merged_coeqtl_df = merged_coeqtl_df.rename({
                'MetaBeta_discovery': 'MetaBeta',
                'MetaSE_discovery': 'MetaSE'
            },
                axis=1)
            merged_coeqtl_df.to_csv(
                workdir / f'output/{filtertype}/rb_calculations/discovery_{celltype_discovery}_replication_{celltype_replication}.fixed.tsv.gz',
                sep='\t',
                compression='gzip')