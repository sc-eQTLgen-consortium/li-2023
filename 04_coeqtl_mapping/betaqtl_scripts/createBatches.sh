condition=$1
celltype=$2
workdir="/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/"
coexpressionfile=${workdir}/"input/individual_networks/${condition}/${condition}_${celltype}.onemillionv23stemiv2ng.zscores.tsv.gz"
gtefile=${workdir}/"input/summary/gte-fix.tsv"
gtfile=${workdir}/"output/genotypevcfs/chrCHR/GenotypeData.bgz.vcf.gz"
batchsize=100000

genelist=${workdir}/"output/${condition}_${celltype}/genelist.noduplicated.txt"
geneannotation=${workdir}/"input/summary/${condition}_${celltype}.genepairs.annotation.gene1position.noduplicated.tsv"
jobtemplatefile=${workdir}/"output/betaqtl_scripts/jobtemplate.noduplicated.sh"
outputfile=${workdir}/"output/${condition}_${celltype}/noduplicated/"
mkdir -p ${outputfile}
python createBatches.py \
	${coexpressionfile} \
	${gtefile} \
	${gtfile} \
	${genelist} \
	${geneannotation} \
	${jobtemplatefile} \
	${batchsize} \
	${outputfile} \
	${condition} \
	${celltype}

genelist=${workdir}/"output/${condition}_${celltype}/genelist.duplicatedversion1.txt"
geneannotation=${workdir}/"input/summary/${condition}_${celltype}.genepairs.annotation.gene1position.duplicatedversion1.tsv"
jobtemplatefile=${workdir}/"output/betaqtl_scripts/jobtemplate.duplicatedversion1.sh"
outputfile=${workdir}/"output/${condition}_${celltype}/duplicatedversion1/"
mkdir -p ${outputfile}
python createBatches.py \
	${coexpressionfile} \
	${gtefile} \
	${gtfile} \
	${genelist} \
	${geneannotation} \
	${jobtemplatefile} \
	${batchsize} \
	${outputfile} \
	${condition} \
	${celltype}

genelist=${workdir}/"output/${condition}_${celltype}/genelist.duplicatedversion2.txt"
geneannotation=${workdir}/"input/summary/${condition}_${celltype}.genepairs.annotation.gene1position.duplicatedversion2.tsv"
jobtemplatefile=${workdir}/"output/betaqtl_scripts/jobtemplate.duplicatedversion2.sh"
outputfile=${workdir}/"output/${condition}_${celltype}/duplicatedversion2"
mkdir -p ${outputfile}
python createBatches.py \
	${coexpressionfile} \
	${gtefile} \
	${gtfile} \
	${genelist} \
	${geneannotation} \
	${jobtemplatefile} \
	${batchsize} \
	${outputfile} \
	${condition} \
	${celltype}
