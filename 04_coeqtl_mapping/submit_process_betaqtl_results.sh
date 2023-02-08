#!/usr/bin/env bash
#SBATCH --time=01:00:00
#SBATCH --mem=20gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module purge

conda init bash
source /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/tools/Beeline/miniconda/etc/profile.d/conda.sh
conda activate scpy3.8


celltype=$1

condition='UT'
workdir="./"

# unfiltered results
python ${workdir}/output/concat_betaqtl_results.fixed.py \
--prefix ${workdir}/output/unfiltered_results/${condition}_${celltype_individual} \
--savepath ${workdir}/output/unfiltered_results/${condition}_${celltype_individual}/concated_alltests_output_fixed.tsv.gz \
--annotation_prefix ${workdir}/input/summary/${condition}_${celltype}
python ${workdir}/output/screen_permutation_p_values.py \
--eqtl_path ${workdir}/input/snp_selection/eqtl/${condition}_${celltype}_eQTLProbesFDR0.05-ProbeLevel.tsv \
--result_prefix ${workdir}/output/unfiltered_results/${condition}_${celltype_individual} \
--save_prefix ${workdir}/output/unfiltered_results/${condition}_${celltype_individual} \
--annotation_prefix ${workdir}/input/summary/${condition}_${celltype}
python ${workdir}/output/multipletesting_correction.fixed.py \
--permutation_pvalue_path ${workdir}/output/unfiltered_results/${condition}_${celltype_individual}/concated_alltests_permutations_fixed.tsv.gz \
--coeqtl_path ${workdir}/output/unfiltered_results/${condition}_${celltype_individual}/concated_alltests_output_fixed.tsv.gz \
--eqtl_path ${workdir}/input/snp_selection/eqtl/${condition}_${celltype}_eQTLProbesFDR0.05-ProbeLevel.tsv \
--save_prefix ${workdir}/output/unfiltered_results/${condition}_${celltype_individual}/coeqtls_fullresults_fixed

# filtered results
python ${workdir}/output/concat_betaqtl_results.fixed.py \
--prefix ${workdir}/output/filtered_results/${condition}_${celltype_individual} \
--savepath ${workdir}/output/filtered_results/${condition}_${celltype_individual}/concated_alltests_output_fixed.tsv.gz \
--annotation_prefix ${workdir}/input/summary/${condition}_${celltype}
python ${workdir}/output/screen_permutation_p_values.py \
--eqtl_path ${workdir}/input/snp_selection/eqtl/${condition}_${celltype}_eQTLProbesFDR0.05-ProbeLevel.tsv \
--result_prefix ${workdir}/output/filtered_results/${condition}_${celltype_individual} \
--save_prefix ${workdir}/output/filtered_results/${condition}_${celltype_individual} \
--annotation_prefix ${workdir}/input/summary/${condition}_${celltype}
python ${workdir}/output/multipletesting_correction.fixed.py \
--permutation_pvalue_path ${workdir}/output/filtered_results/${condition}_${celltype_individual}/concated_alltests_permutations_fixed.tsv.gz \
--coeqtl_path ${workdir}/output/filtered_results/${condition}_${celltype_individual}/concated_alltests_output_fixed.tsv.gz \
--eqtl_path ${workdir}/input/snp_selection/eqtl/${condition}_${celltype}_eQTLProbesFDR0.05-ProbeLevel.tsv \
--save_prefix ${workdir}/output/filtered_results/${condition}_${celltype_individual}/coeqtls_fullresults_fixed
