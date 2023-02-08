"""
RNA velocity analysis using the dynamic model
run on all samples of Oelen v3 dataset for classical monocytes (mon1, mono2)
and filtered for the 2000 most variable genes

Input: loom files generated from velocyto
Output: hd5ad object with RNA velocity estimates
"""

import scvelo as scv
import pandas as pd
import os

scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params('scvelo')  # for beautified visualization

#Load annotation file with UMAP coordinates
fpath="annotations/umap_monocytes.tsv"
umap_coords=pd.read_csv(fpath, sep='\t')

#Load data for each processed lane
lanes=os.listdir("velocyto")
ldata_array=[]
for lane in lanes:
    
    print(lane)
    
    #Get loom file for each lane (file name unfortnuately not always the same)
    files=os.listdir("velocyto/"+lane)
    file=[f for f in files if f.endswith(".loom")]
    lfile="velocyto/"+lane+"/"+file[0]
    
    #Read file
    ldata = scv.read(lfile, cache=True)
    
    #Filter monocytes from file
    filteredNames=[barcodeName.split(":")[1] for barcodeName in ldata.obs.index]
    filteredNames=[barcodeName.replace("x","")+"_"+lane for barcodeName in filteredNames]
    ldata.obs.index=filteredNames
    umap_coords_filtered=umap_coords[umap_coords["Unnamed: 0"].isin(filteredNames)]
    
    #Filter for monocyotes (barcode in umap file)
    ldata=ldata[umap_coords_filtered["Unnamed: 0"],:].copy()
    
    #Make variable names unique
    ldata.var_names_make_unique()
    #Add ldata object
    ldata_array.append(ldata)

ldata_filtered=ldata_array[0].concatenate(ldata_array[1:], batch_key='lane', 
                          batch_categories=lanes,index_unique=None)

#Delete variables which are not required anymore
del ldata
del ldata_array

#Add information about cell types and time points (more annotations are available)
annotations=pd.read_csv("seurat_object_meta.tsv", sep='\t')
annotations=annotations.loc[ldata_filtered.obs.index.tolist(),:]
ldata_filtered.obs["timepoint"]=annotations["timepoint"]
ldata_filtered.obs["celltype"]=annotations["cell_type"]

#Filter for only monocytes 1 and 2
ldata_filtered=ldata_filtered[ldata_filtered.obs.celltype.isin(["mono 1","mono 2"])].copy()

#RNA velocity analysis
scv.utils.show_proportions(ldata_filtered)
#Filter genes with less than 20 counts (spliced + unspliced) and 
#reduce to the top 2000 highly variable genes
scv.pp.filter_and_normalize(ldata_filtered, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(ldata_filtered)

#Run dynamic model
scv.tl.recover_dynamics(ldata_filtered)
scv.tl.velocity(ldata_filtered, mode='dynamical')
scv.tl.velocity_graph(ldata_filtered)

#Add UMAP coordinates
umap_coords.index=umap_coords["Unnamed: 0"]
umap_coords=umap_coords.loc[ldata_filtered.obs.index.tolist(),:]
ldata_filtered.obsm["X_umap"]=umap_coords[["umap_1","umap_2"]].to_numpy()

#Save file
ldata_filtered.write("h5ad_objects/scveloAnalysis_dynamic_velocity_womono34.h5ad")

#Create plot with embedding
scv.pl.velocity_embedding_stream(ldata_filtered, basis='umap', color=['timepoint', 'celltype'],
                                 show=False,save="embedding_dynamic_monocytes_womono34.png")

scv.pl.velocity_graph(ldata_filtered,color="timepoint",
                       show=False,save="velocityGraph_dynamic_monocytes_womono34.png")

#Calculate pseudotime
scv.tl.latent_time(ldata_filtered)
scv.pl.scatter(ldata_filtered, color='latent_time', cmap='gnuplot',
               show=False,save="latenttime_dynamic_monocytes_womono34.png")

#Save file
ldata_filtered.write("h5ad_objects/scveloAnalysis_dynamic_velocity_latenttime_womono34.h5ad")
