# 04_coeqtl_mapping

*plot_effect_concordance_across_cohorts.R*:  compares effect sizes (Z-scores) calculated in each individual dataset (before the meta-analysis)

*plot_celltype_overlap_upset.R* : upset plot overlap of significant co-eQTLs between cell types

*power_analyis_coeqtls.R* : explores how the number of tests reduce the power to detect co-eQTLs, taking estimates for number of tests based on how many genes are expressed above different cutoffs for Oelen v3 dataset

Rb calculations are these files:
*prepare_for_rb_calculation.py* : prepare the input files for rb calculation <br />
*Rb.R* : rb function <br />
*calculate_rb_for_sc_and_bios.R* : execute rb functions<br />
*rb_celltypes.ipynb*: examine the rb values for each cell type; also include scripts for examining different characteristics for coeQTLs compared to non-coeQTLs<br />

Co-eQTL pipeline are these files:<br />
all files in the betaqtl_scripts (incl. templates)<br />
*individual_networks.py*: make co-expression files for each individual<br />
*prepare_genelist_and_annotation_for_betaqtl.py*: prepare input files for the qtl mapping pipeline<br />
*createBatches.sh*: create batches for qtl mapping pipeline<br />
*submit_process_betaqtl_results.sh*: submit the jobs for concatenating qtl mapping and perform multiple testing procedures<br />
*concat_betaqtl_results.fixed.py*: concat qtl mapping results<br />
*screen_permutation_p_values.py*: concat permutation files<br />
*multipletesting_correction.fixed.py*: perform multiple testing correction<br />

Other co-eQTL analysis:<br />
*filtering_strategy.py*: filter for gene pairs<br />
*individual_networks_cmono_ncmono.py*: create co-expression files for each individual for sub cell types in monocytes<br />
*individual_networks_maxcell.py*: create co-expression files for each individual with a limit of cell number<br />
*merge_coexpression_for_betaeqtl_maxcell.py*: merge the co-expression files for each individual with a limit of cell number<br />
*merge_coexpression_for_betaqtl.subsampleindividuals.py*: create co-expression files for each individual with a limit of sample number<br />

BIOS replication are these files:<br />
*replication_in_bios.py*: perform bios replication<br />
*select_snps_from_vcf.sh*: select SNP from vcf file<br />
*examine_bios_replication.ipynb*: examine the bios replication results<br />

Annotating coeQTL results:<br />
*annotate_coeqtl_files.py*: annotate the coeqtl results for nonzero ratio, mean and var of gene pair<br />
*collect_nonzeroratio.py*: collect non zero ratio annotation for all genes in all datasets<br />





