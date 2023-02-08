import pandas as pd
from pathlib import Path

def find_eqtlsnp_gene(snp_genepair, eqtl_snp_gene_set):
    snp = snp_genepair.split('_')[0]
    gene1, gene2 = snp_genepair.split('_')[1].split(';')
    if '_'.join([snp, gene1]) in eqtl_snp_gene_set:
        return gene1
    else:
        return gene2


workdir = Path("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/output")
eqtl_prefix = Path("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/input/snp_selection/eqtl")
writer = pd.ExcelWriter(workdir/'summary/coeQTLs_6majorcelltypes.unfiltered.xlsx', engine='xlsxwriter')
for celltype in ['CD4T', 'CD8T', 'monocyte', 'B', 'DC', 'NK']:
    eqtls_path = eqtl_prefix/f'UT_{celltype}_eQTLProbesFDR0.05-ProbeLevel.tsv'
    eqtl_df = pd.read_csv(eqtls_path, sep='\t')
    eqtl_df['snp_gene'] = ['_'.join(item) for item in eqtl_df[['SNPName', 'genename']].values]
    eqtl_snp_gene_set = set(eqtl_df['snp_gene'])
    df = pd.read_csv(workdir/f'unfiltered_results/UT_{celltype}/coeqtls_fullresults.sig.tsv.gz', sep='\t', compression='gzip')
    df['eqtlgene'] = [find_eqtlsnp_gene(item, eqtl_snp_gene_set) for item in df['snp_genepair']]
    print(celltype, df.shape[0], len(df['eqtlgene'].unique()))
    df.to_excel(writer, sheet_name=celltype)
writer.save()


writer = pd.ExcelWriter(workdir/'summary/coeQTLs_6majorcelltypes.filtered.xlsx', engine='xlsxwriter')
for celltype in ['CD4T', 'CD8T', 'monocyte', 'B', 'DC', 'NK']:
    eqtls_path = eqtl_prefix/f'UT_{celltype}_eQTLProbesFDR0.05-ProbeLevel.tsv'
    eqtl_df = pd.read_csv(eqtls_path, sep='\t')
    eqtl_df['snp_gene'] = ['_'.join(item) for item in eqtl_df[['SNPName', 'genename']].values]
    eqtl_snp_gene_set = set(eqtl_df['snp_gene'])
    df = pd.read_csv(workdir/f'filtered_results/UT_{celltype}/coeqtls_fullresults.sig.tsv.gz', sep='\t', compression='gzip')
    df['eqtlgene'] = [find_eqtlsnp_gene(item, eqtl_snp_gene_set) for item in df['snp_genepair']]
    print(celltype, df.shape[0], len(df['eqtlgene'].unique()))
    df.to_excel(writer, sheet_name=celltype)
writer.save()