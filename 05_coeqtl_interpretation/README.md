# 05_coeqtl_interpretation

*enrichment_GO_terms.R* :  checks GO enrichment among co-eGenes

*enrichment_GO_terms_reduced_background.R* : alternative version of *enrichment_GO_terms.R*  with a more specific background, including only genes tested for the respective SNP-eGene combination

*enrichment_TFs_Remap_preprocessing.R*: filters TFBS annotations from Remap2022 for blood-related cell lines

*enrichment_TFs_Remap.R* : checks enrichment of TFBS among co-eGenes using Remap 2022 annotations in three steps: 1) checks for each cell type and coeGene cluster the enrichment (FDR corrected), 2) checks if the enriched TF is itself part of the coeGenes, 3) checks if the SNP or a SNP in LD is part of the TF

*general_stats.R*: shows general distribution of co-eQTLs (coeGenes per eQTL, distribution of direction of effect)

*magma_coeqtl.Rmd*: R-markdown to perform MAGMA analysis for  GWAS enrichment separately for each set of coeGenes that share the same eQTL (for all eQTLs with at least 5 coeGenes)

*plot_CD4T_mono_network.R*: Create network of all co-eQTLs from CD4+ T cells and/or Monocytes (nodes: eQTLs and coeGenes), color edges by direction of effect

*plot_CD4T_mono_RPS26_subnetwork.R*: Part of the co-eQTL network from *plot_CD4T_mono_network.R* connected with rs1131017-RPS26

*snipe.R*: help functions to query SNiPA website for SNPs in high LD with query SNPs, used in *enrichment_TFs_Remap.R*

*LDTRAIT.ipynb*: Add GWAS annotation to SNP / SNP in LD

*TEM_NAIVE.ipynb*: Examine the impact of rs1131017 on the ratio between TEM and naive T cells

*MS1_Libraries.r*: libraries used for R1_TRANSFAC_enrichment.ipynb, R2_Coloc.ipynb, R3_Coloc_Evaluation.ipynb scripts

*R1_TRANFAC_enrichment.ipynb*: checks TRANSFAC enrichment among co-eGenes and compares them to the Remap enrichment results

*R2_Coloc.ipynb*: runs colocalization analysis on the output of the eqtl and co-eqtl pipeline

*R3_Coloc_Evaluation.ipynb*: investigates and summarizes the results from the colocalization analysis of script R2_Coloc.ipynb


