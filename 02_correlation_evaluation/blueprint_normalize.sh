#!/usr/bin/env bash
#SBATCH --time=7:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=30gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module purge
module load Java

jar_file=eqtl-mapping-pipeline-1.4.9-SNAPSHOT/eqtl-mapping-pipeline.jar
traitfile=./blueprint/tcel_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.txt.gz
outdir=./blueprint
logFile=./blueprint/blueprint_cd4t_adjustPCA.log
java -Xmx30g -Xms30g -jar ${jar_file} \
--mode normalize \
--in ${traitfile} \
--out ${outdir} \
--adjustPCA \
--maxnrpcaremoved 3 \
--stepsizepcaremoval 1 | tee ${logFile}

jar_file=eqtl-mapping-pipeline-1.4.9-SNAPSHOT/eqtl-mapping-pipeline.jar
traitfile=./blueprint/mono_gene_nor_combat_20151109.ProbesWithZeroVarianceRemoved.ProbesCentered.SamplesZTransformed.txt.gz
outdir=./blueprint
logFile=./blueprint/blueprint_normalize.log
java -Xmx30g -Xms30g -jar ${jar_file} \
--mode normalize \
--in ${traitfile} \
--out ${outdir} \
--adjustPCA \
--maxnrpcaremoved 3 \
--stepsizepcaremoval 1 | tee ${logFile}