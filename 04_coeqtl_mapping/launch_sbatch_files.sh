# Calculate individual networks
working_dir=/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping
condition='UT'
for celltype in 'CD4T' 'CD8T' 'B' 'NK' 'DC'
do
for dataset in 'stemiv2' 'onemillionv2' 'onemillionv3' 'ng'
do
  echo ${dataset}_${condition}_${celltype}
  sbatch --parsable --job-name ${dataset}_${condition}_${celltype} \
  --output ${working_dir}/input/individual_networks/logs/${dataset}_${condition}_${celltype}.out \
  --error ${working_dir}/input/individual_networks/logs/${dataset}_${condition}_${celltype}.err \
  ${working_dir}/input/individual_networks/submit_individual_networks.sh ${dataset} ${celltype} ${condition}
done
done # decided not to save into tsv after saving in numpy

# merge individual networks and create gene list and annotation file for betaqtl
for celltype in 'CD4T' 'CD8T' 'B' 'NK' 'DC'
do
  echo ${condition}_${celltype}
  sbatch --parsable --job-name merge_${condition}_${celltype} \
  --output ${working_dir}/input/individual_networks/logs/merge_${condition}_${celltype}.out \
  --error ${working_dir}/input/individual_networks/logs/merge_${condition}_${celltype}.err \
  ${working_dir}/input/individual_networks/submit_merge_coexpression.sh ${celltype} ${condition}
done


# rsync the betaqtl_scripts to gearshift: /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/output/${condition}_${celltype}
# make batches for betaqtl
for celltype in 'CD4T' 'CD8T' 'B' 'NK' 'DC'
do
cd ${working_dir}/output/${condition}_${celltype} ||exit
./createBatches.sh ${condition} ${celltype}
# submit betaqtl jobs
./suball.sh ${working_dir}/output/${condition}_${celltype}/noduplicated/jobs
./suball.sh ${working_dir}/output/${condition}_${celltype}/duplicatedversion1/jobs
./suball.sh ${working_dir}/output/${condition}_${celltype}/duplicatedversion2/jobs
done

# concate and process output from betaqtl
for celltype in 'CD4T' 'CD8T' 'B' 'NK' 'DC'
do
  cd ${working_dir}/output/${condition}_${celltype} ||exit
  echo ${condition}_${celltype}
  sbatch --parsable --job-name process_betaqtl_results_${condition}_${celltype} \
  --output ${working_dir}/input/individual_networks/logs/process_betaqtl_results_${condition}_${celltype}.out \
  --error ${working_dir}/input/individual_networks/logs/process_betaqtl_results_${condition}_${celltype}.err \
  ${working_dir}/output/submit_process_betaqtl_results.sh ${condition} ${celltype}
done