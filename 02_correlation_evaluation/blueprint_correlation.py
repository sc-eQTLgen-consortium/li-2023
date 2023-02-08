import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr, pearsonr
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
# %matplotlib inline

# %%bash
# export HDF5_USE_FILE_LOCKING='FALSE'

def read_numpy(fileprefix, rowname='rows', colname='cols'):
    data = np.load(fileprefix+'.npy')
    rows = [item.strip() for item in open(fileprefix+f'.{rowname}.txt', 'r').readlines()]
    cols = [item.strip() for item in open(fileprefix+f'.{colname}.txt', 'r').readlines()]
    return pd.DataFrame(data=data,
                        index=rows,
                        columns=cols)

def get_pairwise_correlations(corr_df):
    corrmatrix = corr_df.corr()
    triuindices = np.triu_indices(corrmatrix.shape[0], k=1)
    return corrmatrix.values[triuindices]

blueprint_mappings = pd.read_csv('../blueprint/blueprint_mappings.txt',
                                 sep='\t', index_col=0)['Gene name'].T.to_dict()
data = pd.read_csv('../blueprint/mono_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.1PCAsOverSamplesRemoved.txt.gz',
                  sep='\t', index_col=0, compression='gzip')

data.index = [item.split('.')[0] for item in data.index]
data['genename'] = [blueprint_mappings.get(ids) for ids in data.index]
print(data.shape)
data = data.dropna(subset=['genename']).drop_duplicates(subset=['genename'])
data = data.set_index('genename')

data.head()
coefs, ps = spearmanr(data, axis=1)
print(coefs.shape)
np.save('mono_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.1PCAsOverSamplesRemoved.spearmanr.npy',
       coefs)
np.save('mono_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.1PCAsOverSamplesRemoved.spearmanrPvalues.npy',
       ps)
with open('mono_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.1PCAsOverSamplesRemoved.spearmanr.genes.txt',
         'w') as f:
    f.write('\n'.join(data.index.values))