# -----------------------------------------------------------------------------
# Count number of second genes per eQTL (across all cell types)
# And the distribution of effect sizes across cell types
# Input: significnat coeQTL results for all cell types
# Output: printed summary statistics
# -----------------------------------------------------------------------------

library(data.table)
library(dplyr)

coeqtl_dir<-"coeqtl_mapping/output/filtered_results/"

all_counts<-NULL
for(cell_type in c("CD4T","CD8T","monocyte","B","NK","DC")){
  
  # Read coeQTL results
  coeqtls<-fread(paste0(coeqtl_dir,"UT_",
                        cell_type,"/coeqtls_fullresults_fixed_withAF.sig.tsv"))
  
  # Correct MAF
  coeqtls$MetaPZ<-ifelse(coeqtls$AF>0.5,(-1)*coeqtls$MetaPZ,coeqtls$MetaPZ)
  
  count_second_genes<-coeqtls%>%
    group_by(snp_eqtlgene)%>%
    summarise(num_second_genes=n(),
              pos_dir=sum(MetaPZ>0),
              neg_dir=sum(MetaPZ<0))
  
  count_second_genes$ct<-cell_type
  all_counts<-rbind(all_counts,count_second_genes)
}

# Get some general statistics about number of coeGenes per eQTL
# summary(all_counts$num_second_genes)
# sum(all_counts$num_second_genes>4)
# mean(all_counts$num_second_genes>4)
unique_eqtls<-length(unique(all_counts$snp_eqtlgene))
unique_eqtls_second<-length(unique(all_counts$snp_eqtlgene[all_counts$num_second_genes>=5]))
unique_eqtls_second/unique_eqtls

# Get general proportion of positive / negative co-eGenes
count_direction<-all_counts%>%
  group_by(ct)%>%
  summarise(pos_dir=sum(pos_dir),
            neg_dir=sum(neg_dir),
            frac_pos=sum(pos_dir)/sum(num_second_genes))

# Filter again for all genes with at least 5 genes
all_counts_filtered<-all_counts[all_counts$num_second_genes >= 5,]
all_counts_filtered$frac_pos<- all_counts_filtered$pos_dir / all_counts_filtered$num_second_genes

