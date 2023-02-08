# 03_celltype_individual_comparison

*compare_individuals_variance.R* : explores for all genes expressed in at least 50% of the cells the variance across individuals

*correlation_between_celltypes.R* : calculates the Pearson correlation of gene pairwise Spearman correlation for all pairwise combinations of cell types within each dataset (for Oelen v2 and v3 dataset, input from *correlation_celltype.py*), taking only genes expressed in 50% of the cells in both cell types; plots results in heatmap afterwards

*correlation_celltype.py* : calculates Spearman correlation for each genepair expressed in 50% of the cells for Oelen dataset (V2) and (V3), separately per cell type, but combing all individuals; provides so the input csv files for *correlation_between_celltypes.R*

*correlation_correlation_distribution_celltypes_and_individuals.R* : combines two basic overview plots: the correlation distribution in each cell type (input from *correlation_celltype.py*) and the concordance of donor-specific correlation (calculates Pearson correlation of gene pairwise Spearman correlation for each combination of individuals within each cell type)

*correlation_subsampling.py* : calculates per donor correlation for each cell type and different numbers of cells for the sample (randomly subsampling to this number of cells), followed by comparison between donors for within the cell type and the subsampling step, taking genepairs expressed in 50% of the cells, using again Oelen v2 and v3 dataset separately

*fit_logcurve_indiv_subsampling_effect.R* : tkes the results from *correlation_subsampling.py* and fitting logarithmic curves for the  relationship between number of cells and concordance between individuals, one per celltype, to better describe this relationship 

*plot_indiv_subsampling_effect.R* : plots results from *correlation_subsampling.py* to show relationship between number of cells and concordance between individuals
