#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=24g
#SBATCH --cpus-per-task=11
#SBATCH -o LOGPREFIX.log
#SBATCH -e LOGPREFIX.err

set -e
set -u




ml Java/11-LTS
# ml Java/11.0.2

# CHROM, BATCHFILE, OUTPREFIX
# EXP, GTE, GENOTYPE
# CONDITION CELLTYPE
threads=11
java -Xmx17g \
	-Djava.util.concurrent.ForkJoinPool.common.parallelism=$threads \
	-Dmaximum.threads=$threads -Dthread.pool.size=$threads \
	-jar /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/tools/BetaQTL-1.0-SNAPSHOT-jar-with-dependencies.jar \
	-m betaqtl \
	--maf 0.1\
	-a /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/input/summary/CONDITION_CELLTYPE.genepairs.annotation.gene1position.duplicatedversion1.tsv \
	-e /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/input/individual_networks/CONDITION/CONDITION_CELLTYPE.onemillionv23stemiv2ng.zscores.tsv.gz \
	-sgl /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/input/snp_genepair_selection/CONDITION_CELLTYPE.baseline.duplicatedversion1.tsv \
	-gl BATCHFILE \
	-g /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/input/summary/gte-fix.tsv \
	-v GENOTYPE \
	--chr CHROM \
	-o OUTPREFIX \
	--perm 100 \
	--outputall \
	 --snplog \
	 --outputallpermutations
