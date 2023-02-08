# ------------------------------------------------------------------------------
# This scripts runs GO enrichment analysis of significant co-eQTLs results 
# (one analysis per co-eGene group associated with the same eQTL).
# It consists of three functions and requires all cell-type specific eQTL results 
# 'coeqtls_fullresults_fixed.all.tsv.gz' (with all genes tested) and 
# significant results 'coeqtls_fullresults_fixed.sig.tsv.gz' to run it. 
# ------------------------------------------------------------------------------

#################################################################
##                          Libraries                          ##
#################################################################

library(enrichplot)
library(stringr)
library(ggplot2)
library(clusterProfiler)

##################################################################
##          Function for enrichment of subset of genes          ##
##################################################################

# enrich_small function run enrichment analysis of test_genes with background set of genes, 
# store it to the table and also as a dotplot to output pdf file.

enrich_small <- function(test_genes, background_genes, mark_gene, cell_type){

  enrich_out <-enrichGO(gene=c(test_genes),
                        OrgDb='org.Hs.eg.db',
                        keyType="SYMBOL",
                        pvalueCutoff = 0.05,
                        #pAdjustMethod = "none",
                        universe = background_genes,
                        ont="all",
                        minGSSize=5)
  
  enrich_out_df <- data.frame(enrich_out)
  
  # save the data if GO terms found
  if(nrow(enrich_out_df)>0){
    enrich_out_df$id<-mark_gene
    dotplot(enrich_out, showCategory=15) + ggtitle(paste(cell_type, mark_gene))
    ggsave(paste(outdir,  cell_type, "_", mark_gene,  ".pdf", sep=''), width = 10, height = 8)
  }
  return(enrich_out_df)
}

# enrichment_for_zscore_subsets function split the table of co-eQTL results 
# based on the Zscore and run enrichment analysis for these subsets. 
# Also, it creates  hubs of genes and run GO enrichmentfor particular gene hub. 
# All the results are stored to the one datatable, visualized into pdf and saved 
# and separated files. 

enrichment_for_zscore_subsets <- function(tab_all, tab_sign, zscore_id, cell_type){
  
  background_genes  <- unique(c(str_split_fixed(tab_all$Gene, ";", 2)))
  
  #Identify eQTL gene and "second gene"
  test_genes <- str_split_fixed(tab_sign$Gene, ";", 2)
  tab_sign$gene1 <- test_genes[,1]
  tab_sign$gene2 <- test_genes[,2]
  tab_sign$eqtlgene<- str_split_fixed(tab_sign$snp_eqtlgene, "_", 2)[,2]
  tab_sign$secondgene<- with(tab_sign,ifelse(gene1==eqtlgene,gene2,gene1))
  tab_sign$gene1 <- NULL
  tab_sign$gene2 <- NULL
  print(paste(cell_type, dim(tab_sign)[1]))
  
  #Perform enrichment analysis over all coeQTL genes together (and all eQTL / second genes)
  
  all_genes <- enrich_small(test_genes, background_genes, mark_gene=paste(zscore_id, 'all', sep="_"),cell_type)
  gene1_genes <- enrich_small(tab_sign$eqtlgene, background_genes, mark_gene=paste(zscore_id, 'eqtlgene'),cell_type)
  gene2_genes <- enrich_small(tab_sign$secondgene, background_genes, mark_gene=paste(zscore_id, 'secondgene'),cell_type)
  enrichment<-rbind(all_genes,gene1_genes,gene2_genes)
  
  count_hubs <- as.data.frame(table(tab_sign$eqtlgene))
  count_hubs$Var1<-as.character(count_hubs$Var1)
  count_hubs <- count_hubs[order(count_hubs$Freq, decreasing = T),]
  
  #Filter all hub genes to contain more than 5 second genes
  count_hubs<-count_hubs[count_hubs$Freq>5,]
  
  for(hub in seq_len(nrow(count_hubs))){
    hub_gene <- count_hubs$Var1[hub]
    second_genes<-tab_sign$secondgene[tab_sign$eqtlgene == hub_gene] 
    hub_enrich <- enrich_small(second_genes, background_genes, 
                               mark_gene=paste( hub_gene, count_hubs$Freq[hub],zscore_id, sep="_"),cell_type)
    enrichment<-rbind(enrichment,hub_enrich)
  }
  
  #Save cell type in the table
  enrichment$cell_type <- cell_type
  
  #Save results
  name_GO<- paste(outdir, "GO_", zscore_id, cell_type, ".tsv", sep='')
  write.table(enrichment, name_GO, sep='\t', quote=F, row.names = F, col.names = T)
  
}

# enrichment_function load the input files and apply previously described functions 
# to selected genes.

enrichment_function <- function(path, cell_type, outdir){
  
  coeqtls<-fread(paste0(path,"UT_",cell_type,"/coeqtls_fullresults_fixed.all.tsv.gz"))
  
  #Load eqtls on which this is based on
  eqtls<-fread(paste0("coeqtl_mapping/input/snp_selection/eqtl/UT_",
                      cell_type,"_eQTLProbesFDR0.05-ProbeLevel_withAF.tsv"))
  eqtls$SNPpair<-paste0(eqtls$SNPName,"_",eqtls$genename)
  coeqtls<-merge(coeqtls,eqtls[,c("SNPpair","SNPType","AlleleAssessed","OverallZScore",
                                  "AF","alt_allele")],
                 by.x="snp_eqtlgene",by.y="SNPpair")
  
  #Swap the Z score
  coeqtls$MetaPZ<-ifelse(coeqtls$AF>=0.5,coeqtls$MetaPZ*(-1),
                         coeqtls$MetaPZ)
  
  tab_all <- coeqtls
  
  coeqtls_sign<-fread(paste0(path,"UT_",cell_type,"/coeqtls_fullresults_fixed.sig.tsv.gz"))
  
  #Load eqtls on which this is based on
  eqtls_sign<-fread(paste0("coeqtl_mapping/input/snp_selection/eqtl/UT_",
                           cell_type,"_eQTLProbesFDR0.05-ProbeLevel_withAF.tsv"))
  eqtls_sign$SNPpair<-paste0(eqtls_sign$SNPName,"_",eqtls_sign$genename)
  coeqtls_sign<-merge(coeqtls_sign,eqtls_sign[,c("SNPpair","SNPType","AlleleAssessed","OverallZScore",
                                                 "AF","alt_allele")],
                      by.x="snp_eqtlgene",by.y="SNPpair")
  
  #Swap the Z score
  coeqtls_sign$MetaPZ<-ifelse(coeqtls_sign$AF>=0.5,coeqtls_sign$MetaPZ*(-1),
                              coeqtls_sign$MetaPZ)
  tab_sign <- coeqtls_sign
  
  enrichment_for_zscore_subsets(tab_all, tab_sign, zscore_id="All_Zscores",cell_type)
  enrichment_for_zscore_subsets(tab_all[tab_all$MetaPZ >0, ], tab_sign[tab_sign$MetaPZ >0, ], zscore_id="Positive_Zscores",cell_type)
  enrichment_for_zscore_subsets(tab_all[tab_all$MetaPZ < 0, ], tab_sign[tab_sign$MetaPZ <0, ], zscore_id="Negative_Zscores",cell_type)
  
}

##################################################################
##                           Analysis                           ##
##################################################################


for(ct in c("CD8T","CD4T",
            "monocyte","NK","B","DC")){
  enrichment_function(path, cell_type=ct, outdir)
}

# Path where coeQTL results are stored 
path <- '/path/to/coeQTLt/results/'
# Path where results of GO enrich will be stored 
outdir <- '/path/to/outputs/'
# for tests
cell_type ="CD8T"

# 
# run_mono <- enrichment_function(path, cell_type="monocyte", outdir)
# run_cd8t <- enrichment_function(path, cell_type="CD8T", outdir)
# run_cd4t <- enrichment_function(path, cell_type="CD4T", outdir)
# run_b <- enrichment_function(path, cell_type="B", outdir)
# run_dc <- enrichment_function(path, cell_type="DC", outdir)
# run_nk <- enrichment_function(path, cell_type="NK", outdir)
