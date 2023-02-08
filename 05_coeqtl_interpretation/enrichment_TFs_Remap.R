# ------------------------------------------------------------------------------
# Check enrichment of TFBS among co-eGenes using Remap 2022 annotations
# - Check for each cell type and coeGene cluster the enrichment (FDR correction)
# - Check if the enriched TF is itself part of the coeGenes
# - Check if the SNP or a SNP in LD is part of the TF
#
# Input: TFBS from Remap2022 (filtered for blood cell types in script 
#        enrichment_TFs_Remap_preprocessing.R), 
#        gene annotation file,
#        coeQTL results (complete result with all tested genes to 
#        define also the background)
# Output: signficant enrichment results with information about overlap 
#         between TF and co-eQTL SNP
# ------------------------------------------------------------------------------

library(rtracklayer)
library(data.table)

source("snipe.R")

peaks<-import("tfbs_enrichment_remap/remap2022_nr_macs2_hg19_v1_0_blood_related.bed")

# Split name into TF and measured cell lines
ann <- t(matrix(unlist(strsplit(values(peaks)[,"name"], ":", fixed=T)), nrow=2))
colnames(ann) <- c("TF", "conditions")
ann <- as.data.frame(ann,stringsAsFactors=FALSE)

values(peaks)<-ann

#Update seqlevel style to match between genotypes
seqlevelsStyle(peaks)<-"NCBI"
peaks<-keepStandardChromosomes(peaks, pruning.mode="coarse")

#Get number of unique TFs
tfs<-as.character(unique(ann$TF))

################################################################################
# Part 1: check if certain TFs are enriched within the co-eGenes
################################################################################

gene_annot<-import("tfbs_enrichment_remap/genes.gtf")
gene_annot<-gene_annot[gene_annot$type =="gene",]
gene_annot<-keepStandardChromosomes(gene_annot, pruning.mode="coarse")

#Read coeQTL file (all tests)
tf_enrichment_combined<-NULL
total_set_tested_eqtls<-NULL
for(ct in c("CD4T","CD8T","monocyte","NK","DC","B")){
  
  coeqtls<-fread(paste0("coeqtl_mapping/output/filtered_results/UT_",ct,
                        "/coeqtls_fullresults_fixed.all.tsv.gz"))
  coeqtls$gene1<-gsub(";.*","",coeqtls$Gene)
  coeqtls$gene2<-gsub(".*;","",coeqtls$Gene)
  coeqtls$eqtlgene<-gsub(".*_","",coeqtls$snp_eqtlgene)
  coeqtls$gene2<-ifelse(coeqtls$gene1 == coeqtls$eqtlgen, coeqtls$gene2,
                        coeqtls$gene1)
  coeqtls$gene1<-NULL
  tested_genes<-unique(c(coeqtls$eqtlgene,coeqtls$gene2))
  
  #Remove trailing .1 (R issue) of missing genes
  missing_genes<-tested_genes[! tested_genes %in% gene_annot$gene_name]
  tested_genes[! tested_genes %in% gene_annot$gene_name]<-gsub("\\.1$","",missing_genes)
  tested_genes<-unique(tested_genes)
  print(paste("Annotation found for x% of the genes:",
              mean(tested_genes %in% gene_annot$gene_name)))
  
  #Read gene position file and determine TSS
  gene_pos<-gene_annot[gene_annot$gene_name %in% tested_genes,]
  mcols(gene_pos)<-data.frame(gene_name=gene_pos$gene_name)
  
  #Get TSS for the genes
  gene_tss<-promoters(gene_pos,upstream=2000,downstream=2000)
  
  #Check for each TF in ReMap overlap with all gene TSS
  tfbs_ann <- sapply(tfs, function(x) overlapsAny(gene_tss,
                                                      peaks[peaks$TF == x]))
  rownames(tfbs_ann)<-gene_tss$gene_name
  
  #Filter for TFs with at least one found binding
  tfbs_ann<-tfbs_ann[,colSums(tfbs_ann)>0]
  
  #Collapse genes with multiple annotations (hit if at least a hit in one annotation)
  tfbs_ann<- apply(tfbs_ann, 2, tapply, rownames(tfbs_ann),
                     max, na.rm=T)
  #Convert it back into a logical matrix
  tfbs_ann<-matrix(as.logical(tfbs_ann),ncol=ncol(tfbs_ann),
                       dimnames=list(rownames(tfbs_ann),colnames(tfbs_ann)))
  
  #Get significant coeQTLs
  coeqtls_sign<-coeqtls[coeqtls$gene2_isSig,]
  occ_eqtl<-as.data.frame(table(coeqtls_sign$snp_eqtlgene))
  occ_eqtl<-occ_eqtl[occ_eqtl$Freq>=5,]
  occ_eqtl$cell_type<-ct
  total_set_tested_eqtls<-rbind(total_set_tested_eqtls,
                                occ_eqtl)
  
  #Perform Fisher's test for the enrichment
  fisher_all_eqtl<-NULL
  for(eqtl in occ_eqtl$Var1){
    
    #Filter significant gene2s with existing promoter annotation
    sign_gene2<-coeqtls_sign$gene2[coeqtls_sign$snp_eqtlgene==eqtl]
    sign_gene2<-intersect(sign_gene2,rownames(tfbs_ann))
    
    #Iterate over each TF
    fisher_res<-NULL
    for(tf in colnames(tfbs_ann)){
      
      counts<-data.frame(tf_binding=c(sum(tfbs_ann[sign_gene2,tf]),sum(tfbs_ann[,tf])),
                         tf_nonbinding=c(sum(!tfbs_ann[sign_gene2,tf]),sum(!tfbs_ann[,tf])))
      
      res_fisher<-fisher.test(counts,alternative="greater")
      
      fisher_res<-rbind(fisher_res,
                        data.frame(celltype=ct,
                                   eqtl,
                                   tf,
                                   is_coeGene = tf %in% sign_gene2,
                                   fisher_pval=res_fisher$p.value,
                                   tf_coeqtl=counts$tf_binding[1],
                                   notf_coeqtl=counts$tf_nonbinding[1],
                                   tf_background=counts$tf_binding[2],
                                   notf_background=counts$tf_nonbinding[2]))
    }
    
    #Multiple testing correction per eQTL
    fisher_res$fisher_fdr<-p.adjust(fisher_res$fisher_pval,method="BH")
    
    fisher_all_eqtl<-rbind(fisher_all_eqtl,fisher_res)
  }
  
  #Check for CD4T specificallly for RPS26 the positive & negative coeGenes separately
  if(ct=="CD4T"){
    eqtl<-"rs1131017_RPS26"
    
    #Test positive coeGenes (MAF not correctly flipped here)
    sign_gene2<-coeqtls_sign$gene2[coeqtls_sign$snp_eqtlgene==eqtl
                                   & coeqtls_sign$MetaPZ < 0]
    sign_gene2<-intersect(sign_gene2,rownames(tfbs_ann))
    
    #Iterate over each TF
    fisher_res<-NULL
    for(tf in colnames(tfbs_ann)){
      
      counts<-data.frame(tf_binding=c(sum(tfbs_ann[sign_gene2,tf]),sum(tfbs_ann[,tf])),
                         tf_nonbinding=c(sum(!tfbs_ann[sign_gene2,tf]),sum(!tfbs_ann[,tf])))
      
      res_fisher<-fisher.test(counts,alternative="greater")
      
      fisher_res<-rbind(fisher_res,
                        data.frame(celltype=ct,
                                   eqtl=paste0(eqtl,"_positive"),
                                   tf,
                                   is_coeGene = tf %in% sign_gene2,
                                   fisher_pval=res_fisher$p.value,
                                   tf_coeqtl=counts$tf_binding[1],
                                   notf_coeqtl=counts$tf_nonbinding[1],
                                   tf_background=counts$tf_binding[2],
                                   notf_background=counts$tf_nonbinding[2]))
    }
    
    #Multiple testing correction per eQTL
    fisher_res$fisher_fdr<-p.adjust(fisher_res$fisher_pval,method="BH")
    fisher_all_eqtl<-rbind(fisher_all_eqtl,fisher_res)
    
    #Test negative coeGenes (MAF not correctly flipped here)
    sign_gene2<-coeqtls_sign$gene2[coeqtls_sign$snp_eqtlgene==eqtl
                                   & coeqtls_sign$MetaPZ > 0]
    sign_gene2<-intersect(sign_gene2,rownames(tfbs_ann))
    
    #Iterate over each TF
    fisher_res<-NULL
    for(tf in colnames(tfbs_ann)){
      
      counts<-data.frame(tf_binding=c(sum(tfbs_ann[sign_gene2,tf]),sum(tfbs_ann[,tf])),
                         tf_nonbinding=c(sum(!tfbs_ann[sign_gene2,tf]),sum(!tfbs_ann[,tf])))
      
      res_fisher<-fisher.test(counts,alternative="greater")
      
      fisher_res<-rbind(fisher_res,
                        data.frame(celltype=ct,
                                   eqtl=paste0(eqtl,"_negative"),
                                   tf,
                                   is_coeGene = tf %in% sign_gene2,
                                   fisher_pval=res_fisher$p.value,
                                   tf_coeqtl=counts$tf_binding[1],
                                   notf_coeqtl=counts$tf_nonbinding[1],
                                   tf_background=counts$tf_binding[2],
                                   notf_background=counts$tf_nonbinding[2]))
    }
    
    #Multiple testing correction per eQTL
    fisher_res$fisher_fdr<-p.adjust(fisher_res$fisher_pval,method="BH")
    
    fisher_all_eqtl<-rbind(fisher_all_eqtl,fisher_res)
  
  }
  
  tf_enrichment_combined<-rbind(tf_enrichment_combined,
                                fisher_all_eqtl[fisher_all_eqtl$fisher_fdr<0.05,])
  
}

table(tf_enrichment_combined$eqtl,tf_enrichment_combined$celltype)
tf_enrichment_combined[tf_enrichment_combined$is_coeGene,]

################################################################################
# Part 2: check if the SNP (a SNP in high LD) is in the TF peak
################################################################################

tf_enrichment_combined$eqtlsnp<-gsub("_.*","",tf_enrichment_combined$eqtl)

#Get all SNPs in LD with the enriched eQTL SNP
enriched_snps<-unique(tf_enrichment_combined$eqtlsnp)

proxies <- snipa.get.ld.by.snp(enriched_snps,
                               rsquare=0.9,
                               population=c('eur'))

snp_grange<-makeGRangesFromDataFrame(proxies,seqnames.field="CHR",
                                     start.field = "POS2",end.field = "POS2")
snp_grange$snp_name<-proxies$RSID
snp_grange$eqtl_snp<-proxies$QRSID

#Iterate over each SNP to check the overlap with the TF
tf_enrichment_combined$snp_tf_overlap<-FALSE
tf_enrichment_combined$snp_name_overlap<-""
for(i in 1:nrow(tf_enrichment_combined)){
  
  #SNPs in high LD with the eQTL SNP
  snp_subset<-snp_grange[snp_grange$eqtl_snp == tf_enrichment_combined$eqtlsnp[i]]
  
  #All peaks of the respective TF
  peaks_tf<-peaks[peaks$TF == tf_enrichment_combined$tf[i]]
  
  overlap_snp<-overlapsAny(snp_subset,peaks_tf)
  
  tf_enrichment_combined$snp_tf_overlap[i]<-any(overlap_snp)
  tf_enrichment_combined$snp_name_overlap[i]<-paste0(snp_subset$snp_name[overlap_snp],collapse=",")
}

#Save LD information and result table
write.table(tf_enrichment_combined,
            file="tfbs_enrichment_remap/tf_remap_enrichment_results.tsv",
            sep="\t",quote=FALSE,row.names=FALSE)
write.table(proxies,
            file="tfbs_enrichment_remap/ld_proxies_with_position.tsv",
            sep="\t",quote=FALSE,row.names=FALSE)

#Filter for SNPs which overlap with the TF
tf_enrichment_overlap<-tf_enrichment_combined[tf_enrichment_combined$snp_tf_overlap,]

table(tf_enrichment_combined$eqtl,tf_enrichment_combined$celltype)
table(tf_enrichment_overlap$eqtl,tf_enrichment_overlap$celltype)
tfs_coegenes<-tf_enrichment_overlap[tf_enrichment_overlap$is_coeGene,c("celltype","eqtl","tf")]
tfs_coegenes[,c("celltype","eqtl","tf")]


