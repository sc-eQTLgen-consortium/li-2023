#!/usr/bin/env bash
#SBATCH --time=16:00:00
#SBATCH --mem=80gb
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module purge

conda init bash
source /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/tools/Beeline/miniconda/etc/profile.d/conda.sh
conda activate scpy3.8


python /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/input/individual_networks/individual_networks.py \
--datasetname $1 \
--celltype $2 \
--condition $3