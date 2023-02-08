import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import os


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--condition', dest = 'condition')
    parser.add_argument('--celltype', dest='celltype')
    return parser

args = parse().parse_args()
condition, celltype = args.condition , args.celltype

# old code for creating the annotation file..
workdir = Path("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/")
eqtl_annotations_path = workdir/f'input/snp_genepair_selection/annotations/{condition}_{celltype}.baseline.annotatedeQTL.tsv'
savepath = workdir/f'output/{condition}_{celltype}/'
if not os.path.isdir(savepath):
    os.mkdir(savepath)

eqtl_annotations = pd.read_csv(eqtl_annotations_path, sep='\t')

annotation_cols = ['Platform', 'ArrayAddress', 'Symbol', 'Chr', 'ChrStart', 'ChrEnd', 'Probe', 'Seq']
gene_annotation_dict = pd.read_csv('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/probeannotation/singleCell-annotation-stripped.tsv',
                                   sep='\t').set_index('Ensembl').T.to_dict()
genename_ensembl_mapping = pd.read_csv(workdir/'../resources/features_v3_reformated_names.tsv',
                                       sep='\t', names=['Ensembl', 'genename']).set_index('genename')['Ensembl'].T.to_dict()


eqtl_annotations['ArrayAddress'] = eqtl_annotations['genepair_sorted']
eqtl_annotations['Symbol'] = eqtl_annotations['genepair_sorted']
eqtl_annotations['Probe'] = eqtl_annotations['genepair_sorted']
eqtl_annotations['Seq'] = 'NNNNNNN'
getchr = lambda x:gene_annotation_dict.get(x)['Chr'] if x in gene_annotation_dict else np.nan
getchrstart = lambda x:int(gene_annotation_dict.get(x)['ChrStart']) if x in gene_annotation_dict else np.nan
getchrend = lambda x:int(gene_annotation_dict.get(x)['ChrEnd']) if x in gene_annotation_dict else np.nan
eqtl_annotations['Platform'] = 'SingleCell'

eqtl_annotations['eqtlgene'] = [item.split(';')[0] for item in eqtl_annotations['eqtlgene1_gene2']]
eqtl_annotations['eqtlgene_ensembl'] = [genename_ensembl_mapping.get(genename) for genename in eqtl_annotations['eqtlgene']]

eqtl_annotations['Chr'] = [getchr(gene) for gene in eqtl_annotations['eqtlgene_ensembl']]
eqtl_annotations['ChrStart'] = [getchrstart(gene) for gene in eqtl_annotations['eqtlgene_ensembl']]
eqtl_annotations['ChrEnd'] = [getchrend(gene) for gene in eqtl_annotations['eqtlgene_ensembl']]
counts = eqtl_annotations['genepair_sorted'].value_counts()
duplicated_genepairs_set = set(counts[counts>1].index.values)
isdup = lambda x:True if x in duplicated_genepairs_set else False
eqtl_annotations['isdup'] = [isdup(genepair) for genepair in eqtl_annotations['genepair_sorted']]
eqtl_annotations[eqtl_annotations['isdup']==False][['snp', 'genepair_sorted']].to_csv(workdir/f'input/snp_genepair_selection/{condition}_{celltype}.baseline.noduplicated.tsv',
                                                     sep='\t', index=False)
eqtl_annotations[eqtl_annotations['isdup']==False][annotation_cols].to_csv(workdir/f'input/summary/{condition}_{celltype}.genepairs.annotation.gene1position.noduplicated.tsv',
                                                     sep='\t', index=False)
eqtl_annotations[eqtl_annotations['isdup']==False][['genepair_sorted']].to_csv(savepath/'genelist.noduplicated.txt', header=None, index=False)


duplicated = eqtl_annotations[eqtl_annotations['isdup']].drop_duplicates(subset=['genepair_sorted'], keep='first')
duplicated[['snp', 'genepair_sorted']].to_csv(workdir/f'input/snp_genepair_selection/{condition}_{celltype}.baseline.duplicatedversion1.tsv',
                                                     sep='\t', index=False)
duplicated[annotation_cols].to_csv(workdir/f'input/summary/{condition}_{celltype}.genepairs.annotation.gene1position.duplicatedversion1.tsv',
                                                     sep='\t', index=False)
duplicated[['genepair_sorted']].to_csv(savepath/'genelist.duplicatedversion1.txt', header=None, index=False)


duplcated_version2 = eqtl_annotations[eqtl_annotations['isdup']].drop_duplicates(subset=['genepair_sorted'], keep='last')
duplcated_version2[['snp', 'genepair_sorted']].to_csv(workdir/f'input/snp_genepair_selection/{condition}_{celltype}.baseline.duplicatedversion2.tsv',
                                                     sep='\t', index=False)
duplcated_version2[annotation_cols].to_csv(workdir/f'input/summary/{condition}_{celltype}.genepairs.annotation.gene1position.duplicatedversion2.tsv',
                                                     sep='\t', index=False)
duplcated_version2[['genepair_sorted']].to_csv(savepath/'genelist.duplicatedversion2.txt', header=None, index=False)

