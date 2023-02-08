#!/usr/bin/env bash
#SBATCH --job-name=B_1m_v2
#SBATCH --output=/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/EMP_mapping/B/err/B_1m_v2.out
#SBATCH --error=/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/EMP_mapping/B/err/B_1m_v2.err
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=64gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

set -e
ml Java/1.8.0_144

java -jar -Xmx40g -Xms20g -XX:StringTableSize=10000019 -XX:MaxPermSize=512m /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/cis_eqtl_single_cell/EMP_mapping_30_11_2021/eqtl-mapping-pipeline-1.4.9a-SNAPSHOT/eqtl-mapping-pipeline.jar --mode metaqtl --settings /groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/EMP_mapping/B/config/1m_v2.xml
