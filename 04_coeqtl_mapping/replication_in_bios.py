import argparse
import os
import subprocess
from pathlib import Path

import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

workdir = Path("./coeqtl_mapping")
bios_exp_path = 'BIOS_NoRNAPhenoNA_NoSexNA_NoMixups_NoMDSOutlier_20RNAseqAlignemntMetrics/data/gene_read_counts_BIOS_and_LLD_passQC.tsv.SampleSelection.ProbesWithZeroVarianceRemoved.TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.txt.gz'

unique_mappingfile = "./resources/features_v3_reformated_names.tsv"
bios_gt_prefix = Path('./genotypes-hrc-imputed-vcf/')
gte_mapping_path = "./coeqtl_mapping/bios/gte.tsv"


def get_snps_from_vcffile(bashfile_path, vcf_path, snps_path, savepath):
    response = subprocess.run([bashfile_path, vcf_path, snps_path, savepath])
    print(response)
    return None


def get_genes_from_gzipfile(expression_path, gene_path, savepath):
    print("Loading exp dataframe...")
    exp_df = pd.read_csv(expression_path, sep='\t', index_col=0, compression='gzip')
    print("Full exp loaded.")
    genes = pd.read_csv(gene_path, sep='\t')['ensembl']
    common_genes = list(set(genes) & set(exp_df.index.values))
    print(f"Selecting {len(common_genes)} to save. {len(genes) - len(common_genes)} genes not found in BIOS")
    selected_exp_df = exp_df.loc[common_genes]
    genes_dic = pd.read_csv(gene_path, sep='\t').set_index('ensembl')['symbol'].T.to_dict()
    common_genes_names = [genes_dic.get(geneid) for geneid in common_genes]
    selected_exp_df.index = common_genes_names
    selected_exp_df.to_csv(savepath, sep='\t')
    print(f"Selected {selected_exp_df.shape[0]} genes in {savepath}.")
    return selected_exp_df


def make_snps_genes_files_for_coeqtls(coeqtl_path):
    significant_coeqtls = pd.read_csv(coeqtl_path, sep='\t', compression='gzip', index_col=0)
    mappings = pd.read_csv(unique_mappingfile, sep='\t', names=['geneid', 'genename', 'type']).set_index('genename')[
        'geneid'].T.to_dict()
    snps = significant_coeqtls['SNP'].unique()
    genes = [ele for item in significant_coeqtls['Gene'].values for ele in item.split(';')]
    genes = list(set([item for item in genes if item]))
    genes_df = pd.DataFrame(data=[[item, mappings.get(item)] for item in genes],
                            columns=['symbol', 'ensembl']).dropna(subset=['ensembl'])
    print(f"Writing {len(snps)} snps and {len(genes)} genes from coeQTLs.")
    with open(f"{str(coeqtl_path)[:-len('.tsv.gz')]}.snps.txt", 'w') as f:
        f.write('\n'.join(snps))
    genes_df.to_csv(f"{str(coeqtl_path)[:-len('.tsv.gz')]}.genes.tsv",
                    sep='\t', index=False)
    return snps, genes_df


def replicate(annotated_coeqtl_path,
              bios_gt_path,
              bios_gene_path,
              saveprefix,
              gte_mapping_path,
              vcf_header_rows=6):
    import warnings
    def find_gene2(eqtlgene, genepair):
        gene1, gene2 = genepair.split(';')
        if eqtlgene == gene1:
            return gene2
        else:
            return gene1
    warnings.simplefilter(action='ignore', category=FutureWarning)
    gte_mapping = pd.read_csv(gte_mapping_path, sep='\t').set_index('gt')['exp'].T.to_dict()
    # transform the GT columns into expression ids
    gt = pd.read_csv(bios_gt_path, skiprows=vcf_header_rows, sep='\t')
    sc_individuals = pd.read_csv(
        './coeqtl_mapping/input/summary/gte-fix.tsv',
        sep='\t'
    )['genotypesampleID']
    # remove LLD individuals
    remove_individuals = list(set(sc_individuals) & set(gt.columns))
    gt = gt.drop(remove_individuals, axis=1)
    gt_snp_set = set(gt['ID'].values)
    # map genotype and expression individual names
    find_name = lambda x: gte_mapping.get(x) if x in gte_mapping else x
    gt = gt.rename({item: find_name(item) for item in gt.columns}, axis=1)
    # load expression data
    exp = pd.read_csv(bios_gene_path, index_col=0, sep='\t', compression='gzip')
    genename_mapping = pd.read_csv(unique_mappingfile, sep='\t', names=['gene_id', 'gene_name']).set_index('gene_id')[
        'gene_name'].T.to_dict()
    exp['genename'] = [genename_mapping.get(geneid) for geneid in exp.index]
    exp = exp.dropna(subset=['genename']).set_index('genename')
    expression_gene_name_set = set(exp.index)
    common_indidvidauls = list(set(exp.columns) & set(gt.columns))
    gt_df = gt.set_index('ID')
    exp_common = exp[common_indidvidauls]
    coeqtls = pd.read_csv(annotated_coeqtl_path, sep='\t', compression='gzip', index_col=0)
    coeqtls['gene1'] = [item.split('_')[1] for item in coeqtls['snp_eqtlgene']]
    coeqtls['gene2'] = [find_gene2(gene1, genepair) for (gene1, genepair) in coeqtls[['gene1', 'Gene']].values]
    coeqtl_pairs = coeqtls[['SNP', 'gene1', 'gene2']].values
    # start solving the interaction models
    i = 0
    res_df = pd.DataFrame()
    for snp, gene1, gene2 in tqdm(coeqtl_pairs):
        if snp in gt_snp_set and gene1 in expression_gene_name_set and gene2 in expression_gene_name_set:  # todo: ESNG id to genename
            i += 1
            gt_selected = gt_df.loc[snp]
            gene1_selected = exp_common.loc[gene1]
            gene2_selected = exp_common.loc[gene2]
            x_df = pd.concat([gt_selected[common_indidvidauls], gene2_selected], axis=1)
            x_df[f'{snp}_dosage'] = [float(item.split(':')[1]) for item in x_df[snp]]
            x_df[f'{snp}_{gene2}'] = x_df[f'{snp}_dosage'] * x_df[gene2]
            X = sm.add_constant(x_df[[f'{snp}_dosage', gene2, f'{snp}_{gene2}']])
            model = sm.OLS(gene1_selected.T, X)
            results_data = model.fit().summary().tables[1].data
            results = pd.DataFrame(data=results_data[1:], columns=results_data[0]).set_index('')
            results['gene1'] = gene1
            results['gene2'] = gene2
            results['assessed_allele'] = gt_selected['ALT']
            results['num_individuals'] = len(common_indidvidauls)
            res_df = pd.concat([res_df, results], axis=0)
            if len(coeqtl_pairs) > 10000:
                if i % 10000 == 0 and i > 1:
                    res_df.to_csv(f"{saveprefix}.part{int(i / 10000)}.tsv", sep='\t')
                    res_df = pd.DataFrame()
                    print(f"results part {int(i / 10000)} has been saved in {saveprefix}.part{int(i / 10000)}.tsv")
    part_ind = 1 + int(i / 10000)
    res_df.to_csv(f"{saveprefix}.part{part_ind}.tsv.gz", sep='\t', compression='gzip')
    print(f"results part {part_ind} has been saved in {saveprefix}.part{part_ind}.tsv")
    return res_df


def make_gte_mapping_file():
    prefix = Path("./tmp03boxy/input/")
    gtm = pd.DataFrame()
    for filename in os.listdir(prefix / "hrcGTM"):
        sub_gtm = pd.read_csv(prefix / f"hrcGTM/{filename}", compression='gzip', sep='\t', names=['gt', 'met'])
        gtm = pd.concat([gtm, sub_gtm], axis=0)
    gtm = gtm.set_index('met')
    mte_path = prefix / 'hrcMTE/CODAM_LLDeep_LLS660Q_LLSOmni_NTR_RS_MTE.txt'
    mte = pd.read_csv(mte_path, sep='\t', names=['met', 'exp']).set_index('met')
    all_mapping = pd.concat([gtm, mte], axis=1)
    gte = all_mapping[['gt', 'exp']]
    gte = gte.dropna()
    gte.to_csv('./coeqtl_mapping/bios/gte.tsv',
               sep='\t', index=False)
    return gte


def examine_replicated_in_bios(replication_res_path, savepath):
    from statsmodels.stats.multitest import multipletests
    bios_replication = pd.read_csv(replication_res_path, sep='\t')
    bios_replication['snp_genepair'] = ['_'.join(item) for item in bios_replication[['Unnamed: 0', 'gene1']].values]
    tobesave = lambda x: True if 'dosage' not in x and 'const' not in x and x.startswith('rs') else False
    bios_replication['isinteractionterm'] = [tobesave(item) for item in bios_replication['snp_genepair']]
    bios_interactions_df = bios_replication[bios_replication['isinteractionterm']]
    bios_interactions_df['corrected_p'] = multipletests(bios_interactions_df['P>|t|'], method='fdr_bh')[1]
    bios_interactions_df.to_csv(savepath, sep='\t')
    significant_res = bios_interactions_df[bios_interactions_df['corrected_p'] <= 0.05]
    print("Significantly replicated coeQTLs: ", significant_res.shape[0])
    return bios_interactions_df


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--saveprefix', type=str, dest='saveprefix')
    parser.add_argument('--coeqtlpath', type=str, dest='coeqtlpath')
    parser.add_argument('--selection', type=str, dest='selection')
    parser.add_argument('--chromosome', type=str, dest='chr')
    parser.add_argument('--replicate', type=str, dest='replicate')
    parser.add_argument('--bios_selected_vcf', type=str, dest='bios_selected_vcf')
    return parser


if __name__ == '__main__':
    # _ = make_gte_mapping_file() # make the GTE files for BIOS
    args = arguments().parse_args()
    print("Arguments:")
    print(args)
    # get snps and genes from coeQTL results
    coeqtl_path = args.coeqtlpath
    savedirectory = Path(args.saveprefix)
    if not os.path.isdir(savedirectory):
        os.makedirs(savedirectory)
    if args.selection == 'snp':
        if not os.path.exists(f"{coeqtl_path[:-len('.tsv.gz')]}.snps.txt"):
            _ = make_snps_genes_files_for_coeqtls(coeqtl_path)
        # get snps from vcf
        chromosome = args.chr
        snp_bashfile_path = workdir / 'bios/select_snps_from_vcf.sh'
        dataname = f"chr{chromosome}"
        bios_gt_path = bios_gt_prefix / f"chr{chromosome}/GenotypeData.vcf.gz"
        snp_savepath = savedirectory / f'bios_selected.chr{chromosome}.vcf'
        get_snps_from_vcffile(snp_bashfile_path, bios_gt_path, f"{coeqtl_path[:-len('.tsv.gz')]}.snps.txt",
                              snp_savepath)
    elif args.selection == 'gene':
        if not os.path.exists(f"{coeqtl_path[:-len('.tsv.gz')]}.snps.txt"):
            _ = make_snps_genes_files_for_coeqtls(coeqtl_path)
        # get genes from gzip files
        gene_savepath = savedirectory / f'bios_selected_gene_expression.tsv'
        _ = get_genes_from_gzipfile(bios_exp_path, f"{coeqtl_path[:-len('.tsv.gz')]}.genes.tsv", gene_savepath)
    if args.replicate:
        # perform replication in bios
        work_prefix = Path("./coeqtl_mapping/")
        bios_gt_path = args.bios_selected_vcf
        vcf_header_rows = 6
        bios_gene_path = bios_exp_path
        saveprefix = savedirectory / 'bios_replication_results.eqtlgene1_gene2'
        replicate(annotated_coeqtl_path=coeqtl_path,
                  bios_gt_path=bios_gt_path,
                  bios_gene_path=bios_gene_path,
                  saveprefix=saveprefix,
                  gte_mapping_path=gte_mapping_path,
                  vcf_header_rows=vcf_header_rows)
        # concatenate replication results saved in parts
        res_df = pd.DataFrame()
        for filename in os.listdir(savedirectory):
            if filename.startswith('bios_replication_results.eqtlgene1_gene2.') and 'part' in filename and filename.endswith('gz'):
                print(filename)
                df = pd.read_csv(savedirectory / filename, sep='\t', compression='gzip')
                res_df = pd.concat([res_df, df], axis=0)
        tobesave = lambda x: True if 'dosage' not in x and 'const' not in x and x.startswith('rs') else False
        res_df['isinteractionterm'] = [tobesave(item) for item in res_df['Unnamed: 0']]
        bios_interactions_df = res_df[res_df['isinteractionterm']]
        bios_interactions_df['snp_genepair'] = ['_'.join([item[0].split('_')[0],
                                                          ';'.join(sorted([item[0].split('_')[1], item[1]]))]) for item
                                                in
                                                bios_interactions_df[['Unnamed: 0', 'gene1']].values]
        bios_interactions_df['corrected_p'] = multipletests(bios_interactions_df['P>|t|'], method='fdr_bh')[1]
        bios_interactions_df.to_csv(savedirectory / 'bios_replication_results.eqtlgene1_gene2.all.tsv.gz',
                                    sep='\t', index=False, compression='gzip')
        bios_interactions_df[bios_interactions_df['corrected_p'] <= 0.05].to_csv(
            savedirectory / 'bios_replication_results.eqtlgene1_gene2.sig.tsv.gz',
            sep='\t', index=False, compression='gzip')
