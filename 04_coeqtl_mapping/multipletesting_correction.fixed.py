import pandas as pd
from statsmodels.stats.multitest import multipletests
import numpy as np
import argparse
from scipy.optimize import minimize
from scipy.stats import beta
from scipy import special
from pathlib import Path


def read_numpy(prefix):
    data = np.load(f'{prefix}.npy')
    columns = [f'perm{item.strip()}' for item in open(f'{prefix}.cols.txt', 'r').readlines()]
    rows = [item.strip() for item in open(f'{prefix}.rows.txt', 'r').readlines()]
    return pd.DataFrame(data=data, columns=columns, index=rows)


def beta_distribution_mle_function(x, p):
    k, n = x
    ll = (k - 1) * np.sum(np.log(p)) + (n - 1) * np.sum(np.log(1 - p)) - np.size(p) * special.betaln(k, n)
    return -1 * ll


def beta_distribution_initial_guess(x):
    """
    https://stats.stackexchange.com/questions/13245/which-is-a-good-tool-to-compute-parameters-for-a-beta-distribution
    """
    mean = np.mean(x)
    var = np.var(x)
    a = mean * ((mean * (1 - mean) / var) - 1)
    b = (1 - mean) * ((mean * (1 - mean) / var) - 1)
    return a, b


def fit_beta_distribution(p, a_bnd=(0.1, 10), b_bnd=(1, 1000000)):
    a, b = beta_distribution_initial_guess(p)
    x0 = np.array([min(max(a, a_bnd[0]), a_bnd[1]), min(max(b, b_bnd[0]), b_bnd[1])])
    res = minimize(beta_distribution_mle_function,
                   x0=x0,
                   args=(p, ),
                   method='nelder-mead',
                   bounds=(a_bnd, b_bnd),
                   options={"maxiter": 10000, "disp": True})
    return res.x, res.nfev, res.nit


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--permutation_pvalue_path', dest='permutation_pvalue_path')
    parser.add_argument('--coeqtl_path', dest='coeqtl_path')
    parser.add_argument('--eqtl_path', dest='eqtl_path')
    parser.add_argument('--save_prefix', dest='saveprefix')
    return parser


def find_eqtlsnp_gene(snp, genepair, eqtl_snp_gene_set):
    gene1, gene2 = genepair.split(';')
    if '_'.join([snp, gene1]) in eqtl_snp_gene_set:
        return '_'.join([snp, gene1])
    else:
        return '_'.join([snp, gene2])


def find_eqtl_gene(coeqtl_chrpos, annotation_dict):
    annotation_eqtlgene = annotation_dict.get(coeqtl_chrpos)
    return annotation_eqtlgene

def main():
    args = arguments().parse_args()
    coeqtl_path = args.coeqtl_path
    eqtls_path = args.eqtl_path
    saveprefix = args.saveprefix
    permutation_pvalue_path = args.permutation_pvalue_path
    permutation_cols = [f'Perm{ind}' for ind in range(0, 100)]
    permutation_pvalues_df = pd.read_csv(permutation_pvalue_path, sep='\t',
                                         compression='gzip', index_col=0)
    eqtl_df = pd.read_csv(eqtls_path, sep='\t')
    eqtl_df['chr_pos'] = ['_'.join([str(ele) for ele in item]) for item in eqtl_df[['ProbeChr', 'ProbeCenterChrPos']].values]
    eqtl_snp_gene_set = set(['_'.join(item) for item in eqtl_df[['SNPName', 'genename']].values])
    annotation_path = '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/probeannotation/singleCell-annotation-stripped.tsv'
    mappingdic = pd.read_csv('/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/resources/features_v3_reformated_names.tsv',
                              sep='\t', names=['geneid', 'genename']).set_index('geneid')['genename'].T.to_dict()
    annotation_df = pd.read_csv(annotation_path, sep='\t')
    annotation_df['chr_pos'] = ['_'.join([str(ele) for ele in item]) for item in annotation_df[['Chr', 'ChrStart', 'ChrEnd']].values]
    annotation_df['genename'] = [mappingdic.get(ensemblid) for ensemblid in annotation_df['Ensembl']]
    annotation_dict = annotation_df.set_index('chr_pos')['genename'].T.to_dict()
    # eqtl_df['snp_gene'] = ['_'.join(item) for item in eqtl_df[['SNPName', 'genename']].values]
    # eqtl_snp_gene_set = set(eqtl_df['snp_gene'])
    coeqtls = pd.read_csv(coeqtl_path, sep='\t', index_col=0, compression='gzip')
    coeqtls['eqtlgene'] = [find_eqtl_gene(chr_pos, annotation_dict) for (chr_pos) in coeqtls['chr_pos']]
    coeqtls['snp_eqtlgene'] = ['_'.join(item) for item in coeqtls[['SNP', 'eqtlgene']].values]
    coeqtls_lowest_nominalP = coeqtls.sort_values(by='MetaP', ascending=True).drop_duplicates(subset=['snp_eqtlgene'])
    coeqtls_lowest_nominalP_dict = coeqtls_lowest_nominalP.set_index('snp_eqtlgene')['MetaP'].T.to_dict()
    permutation_pvalues_df['SNP'] = [item.split('_')[0] for item in permutation_pvalues_df.index]
    permutation_pvalues_df['nominalP'] = [coeqtls_lowest_nominalP_dict.get(snp) for snp in
                                          permutation_pvalues_df.index]
    permutation_pvalues_df = permutation_pvalues_df.dropna(subset=['nominalP'])
    permutation_pvalues_df['beta_shape1'], permutation_pvalues_df['beta_shape2'] = \
        zip(*[fit_beta_distribution(x)[0] for x in permutation_pvalues_df[permutation_cols].values])
    permutation_pvalues_df['pval_beta'] = [1-beta.sf(x[0], x[1], x[2]) for x in
                                               permutation_pvalues_df[['nominalP', 'beta_shape1', 'beta_shape2']].values]
    assert permutation_pvalues_df['pval_beta'].isnull().sum() == 0
    # over all eqtls, perform BH-FDR
    permutation_pvalues_df['qval'] = multipletests(permutation_pvalues_df['pval_beta'].values, method='fdr_bh')[1]
    permutation_pvalues_df.to_csv(f'{saveprefix}.eqtls_betaadjustedPs.tsv.gz', sep='\t', compression='gzip')
    ub = permutation_pvalues_df[permutation_pvalues_df['qval']>=0.05].sort_values(by=['pval_beta'], ascending=True)['pval_beta'].values[0]
    lb = permutation_pvalues_df[permutation_pvalues_df['qval']<=0.05].sort_values(by=['pval_beta'], ascending=False)['pval_beta'].values[0]
    pthreshold = (ub + lb) / 2
    print('Minimum p-value threshold', pthreshold)
    permutation_pvalues_df['threshold_per_betadistribution'] = [beta.ppf(pthreshold, x[0], x[1]) for x in
                                                                permutation_pvalues_df[['beta_shape1', 'beta_shape2']].values]
    permutation_pvalue_threshold_dict = permutation_pvalues_df.T.to_dict()
    coeqtls['snp_beta_shape1'] = [permutation_pvalue_threshold_dict.get(snp)['beta_shape1'] for snp in coeqtls['snp_eqtlgene'].values]
    coeqtls['snp_beta_shape2'] = [permutation_pvalue_threshold_dict.get(snp)['beta_shape2'] for snp in coeqtls['snp_eqtlgene']]
    coeqtls['snp_pvalbeta'] = [permutation_pvalue_threshold_dict.get(snp)['pval_beta'] for snp in coeqtls['snp_eqtlgene']]
    coeqtls['snp_qval'] = [permutation_pvalue_threshold_dict.get(snp)['qval'] for snp in coeqtls['snp_eqtlgene']]
    coeqtls['gene2_pthreshold'] = [permutation_pvalue_threshold_dict.get(snp)['threshold_per_betadistribution']
                                   for snp in coeqtls['snp_eqtlgene']]
    issig = lambda x:True if x[0] <= x[1] else False
    coeqtls['gene2_isSig'] = [issig(item) for item in coeqtls[['MetaP', 'gene2_pthreshold']].values]
    significant_coeqtls = coeqtls[(coeqtls['snp_qval']<=0.05) & (coeqtls['gene2_isSig'])]
    print('Significant results:', significant_coeqtls.shape[0])
    coeqtls.to_csv(f'{saveprefix}.all.tsv.gz', sep='\t', compression='gzip')
    significant_coeqtls.to_csv(f'{saveprefix}.sig.tsv.gz', sep='\t', compression='gzip')
    return coeqtls


if __name__ == '__main__':
    _ = main()
