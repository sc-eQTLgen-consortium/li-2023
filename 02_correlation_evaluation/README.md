# 02_correlation_evaluation

*blueprint_normalize.sh* normalize BLUEPRINT dataset, as well as regress out the first PC

*blueprint_correlation.py*: calculate the co-expression for gene pairs in BLUEPRINT data

*compare_blueprint_cutoffs_CD4T.py* : Compare correlation between Blueprint and single cell (Oelen v3 dataset) for different expression thresholds (number of cells expressing the gene), implemented for UT and CD4+ T cells here

*compare_immunexut_cutoffs_CD4T.py*: Same approach as in *compare_blueprint_cutoffs_CD4T.py*, but comparing correlation between ImmuNexUT and Oelen v3 dataset instead

*correlation_between_datasets.R*: check Pearson correlation between data sets (for CD4+ T cells) for single cell vs single cell dataset comparison, single cell vs bulk dataset comparison and bulk vs bulk dataset comparison, afterwards combines all results in one large heatmap

*correlation_between_datasets_extended.R*: check if the correlation values between matched cell types for single cell and bulk (ImmuNexUT) are higher than for not-matched cell types

*correlation_between_datasets_othercts.R*: extension of *correlation_between_datasets.R* that includes all cell types (not only CD4+ T cells)

*correlation_timepoint_combined_indivs_1mio.py*: calculate the co-expression for genes that are expressed in more than 50% cells in Oelen v2 and v3 dataset

*correlation_timepoint_combined_indivs_ng.py*: calculate the co-expression for genes that are expressed in more than 50% cells in van der Wijst dataset

*correlation_timepoint_combined_indivs_stemiv2.py*: calculate the co-expression for genes that are expressed in more than 50% cells in van Blokland v2 dataset

*correlation_timepoint_combined_indivs_stemiv3.py*: calculate the co-expression for genes that are expressed in more than 50% cells in van Blokland v3 dataset

*figure2_barplot_cutoffs.R*: create barplots from the results of *compare_blueprint_cutoffs_CD4T.py* and *compare_immunexut_cutoffs_CD4T.py*

*figure2_scatterplots.R*: creates inset plots for Main Figure 2 (a,b,d), showing scatterplots of gene pair-wise Spearman correlation values between two data sets for a) Oelen v3 dataset vs van Blokland v2 dataset (both CD4+ T cells), b) ImmuNexUT - van Blokland v2 (naive CD4+ T cells and CD4+ T cells) and c) Blueprint - ImmuNexUT (both naive CD4+ T cells)

*normalize_ImmuNexUT.R*: preprocessingImmuNexUT data (separately for each cell type with a matching single-cell cell type) following the description in the corresponding publication (filtering lowly expressed genes, TMM normalization and batch correction) followed by correlation calculation for all genes expressed in 50% of the cells of the Oelen v3 dataset (for comparison with single cell data)

*wilcoxon_test_crispr.R*: Benchmark our correlation results from single cell (Oelen v3, CD4+ T cells) and bulk (ImmuNexUT, naive CD4+ T cells) with a public CRISPR perturbation dataset using Wilcoxon Rank Sum Test

*wilcoxon_test_string.R*: Compare if correlated pairs from single cell (Oelen v3, CD4+ T cells) and bulk (ImmuNexUT, naive CD4+ T cells) are enriched in STRING database (Using the same strategy as in CRISPR validation with Wilcoxon Rank Sum Test)


