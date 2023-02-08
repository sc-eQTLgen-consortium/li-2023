#!/usr/bin/env bash
#SBATCH --time=8:00:00
#SBATCH --mem=80gb
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module purge

conda init bash
source /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/tools/Beeline/miniconda/etc/profile.d/conda.sh
conda activate scpy3.8


python /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/input/individual_networks/merge_coexpression_for_betaeqtl.py \
--celltype $1 \
--condition $2

python /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/input/individual_networks/prepare_genelist_and_annotation_for_betaqtl.py \
--celltype $1 \
--condition $2