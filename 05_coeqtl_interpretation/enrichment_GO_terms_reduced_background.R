# ------------------------------------------------------------------------------
# Check for each co-eQTL with at least 5 co-eGenes if there is any enrichment
# using all genes correlated with the respective eGene as background
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(clusterProfiler)

path<-"coeqtl_mapping/output/filtered_results/"
outdir<-"coeqtl_interpretation/go_enrichment/"

#Run the GO enrichment for each cell type
enrichment<-NULL
enrichment_summary<-NULL
coegenes_counts_total<-NULL
for(cell_type in c("CD4T","CD8T","monocyte","NK","B","DC")){

  coeqtls <- fread(paste0(path, "UT_",cell_type, 
                         "/coeqtls_fullresults_fixed.all.tsv.gz"))
  coeqtls$gene1<-gsub(";.*","",coeqtls$Gene)
  coeqtls$gene2<-gsub(".*;","",coeqtls$Gene)
  coeqtls$second_gene<-ifelse(coeqtls$gene1 == coeqtls$eqtlgen, coeqtls$gene2,
                        coeqtls$gene1)
  coeqtls$gene1<-NULL
  coeqtls$gene2<-NULL
  
  # Take all tested genes as background
  background_genes  <- union(coeqtls$eqtlgen,coeqtls$second_gene)
  coeqtls_sign<-coeqtls[coeqtls$gene2_isSig,]
  
  print(paste(cell_type,"with",nrow(coeqtls_sign),"co-eQTLs"))
  print(paste("Size of the combined background set:",
              length(background_genes)))
  
  # Identify all eQTLs with at least 5 coeGenes
  coegene_count<-coeqtls_sign%>%
    group_by(snp_eqtlgene)%>%
    summarize(count_coeGenes=n())%>%
    filter(count_coeGenes>4)
  
  coegene_count$cell_type<-cell_type
  coegenes_counts_total<-rbind(coegenes_counts_total,
                               coegene_count)
  
  #Size of the reduced background set
  coeqtls_reduced_background<-coeqtls%>%
    filter(snp_eqtlgene %in% coegene_count$snp_eqtlgene)%>%
    group_by(snp_eqtlgene)%>%
    summarize(count_secondgene=n())
  
  print("Size of the reduced gene sets")
  print(summary(coeqtls_reduced_background$count_secondgene))

  enrichment_found<-0
  #Perform GO enrichemt separately for each eQTL
  for(eqtl in coegene_count$snp_eqtlgene){
    
    # Run enrichment analysis with background set
    enrich_out <-enrichGO(gene=coeqtls_sign$second_gene[coeqtls_sign$snp_eqtlgene == eqtl],
                          OrgDb='org.Hs.eg.db',
                          keyType="SYMBOL",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          universe = coeqtls$second_gene[coeqtls$snp_eqtlgene == eqtl],
                          ont="all",
                          minGSSize=5)
    
    if(nrow(enrich_out@result)>0){
      
      # Save if a enrichment was found
      enrichment_found<-enrichment_found+1
      
      # Save result dataframe
      res<-enrich_out@result
      res$cell_type<-cell_type
      res$snp_eGene<-eqtl
      enrichment<-rbind(enrichment,
                        res[,c("cell_type","snp_eGene","ONTOLOGY","ID",
                               "Description","pvalue","p.adjust","GeneRatio","BgRatio")])
    }

  }
  
  enrichment_summary<-rbind(enrichment_summary,
                            data.frame(cell_type,
                                       n_eqtls_freq=nrow(coegene_count),
                                       n_enrich=enrichment_found,
                                       freq_enrich=enrichment_found/nrow(coegene_count)))
  
  
  #Check for CD4T specificallly for RPS26 the positive & negative coeGenes separately
  if(cell_type=="CD4T"){
    eqtl<-"rs1131017_RPS26"
    
    #Test positive coeGenes (MAF not correctly flipped here)
    enrich_out <-enrichGO(gene=coeqtls_sign$second_gene[coeqtls_sign$snp_eqtlgene == eqtl &
                                                          coeqtls_sign$MetaPZ < 0],
                          OrgDb='org.Hs.eg.db',
                          keyType="SYMBOL",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          universe = coeqtls$second_gene[coeqtls$snp_eqtlgene == eqtl],
                          ont="all",
                          minGSSize=5)
    
    if(nrow(enrich_out@result)>0){
      
      # Save if a enrichment was found
      enrichment_found<-enrichment_found+1
      
      # Save result dataframe
      res<-enrich_out@result
      res$cell_type<-cell_type
      res$snp_eGene<-paste0(eqtl,"_positive")
      enrichment<-rbind(enrichment,
                        res[,c("cell_type","snp_eGene","ONTOLOGY","ID",
                               "Description","pvalue","p.adjust","GeneRatio","BgRatio")])
    }
    
    #Test negative coeGenes (MAF not correctly flipped here)
    enrich_out <-enrichGO(gene=coeqtls_sign$second_gene[coeqtls_sign$snp_eqtlgene == eqtl &
                                                          coeqtls_sign$MetaPZ > 0],
                          OrgDb='org.Hs.eg.db',
                          keyType="SYMBOL",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          universe = coeqtls$second_gene[coeqtls$snp_eqtlgene == eqtl],
                          ont="all",
                          minGSSize=5)
    
    if(nrow(enrich_out@result)>0){
      
      # Save if a enrichment was found
      enrichment_found<-enrichment_found+1
      
      # Save result dataframe
      res<-enrich_out@result
      res$cell_type<-cell_type
      res$snp_eGene<-paste0(eqtl,"_negative")
      enrichment<-rbind(enrichment,
                        res[,c("cell_type","snp_eGene","ONTOLOGY","ID",
                               "Description","pvalue","p.adjust","GeneRatio","BgRatio")])
    }
  }
}


#Format p-values
enrichment$pvalue<-format(enrichment$pvalue,digits=3)
enrichment$p.adjust<-format(enrichment$p.adjust,digits=3)
write.table(enrichment,
            file=paste0(outdir,"GOenrichment_coeGenes_allcelltypes_otherbackground.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)

#Check general statistics (per eQTL - cell type)
sum(enrichment_summary$n_eqtls_freq)
sum(enrichment_summary$n_enrich)

#Check general statistics (per eQTL, combining all cell types)
enriched_eqtls<-setdiff(unique(enrichment$snp_eGene),
                        c("rs1131017_RPS26_positive","rs1131017_RPS26_negative"))
length(enriched_eqtls)

coegenes_counts_total<-as.data.frame(coegenes_counts_total)
length(unique(coegenes_counts_total$snp_eqtlgene))
setdiff(unique(coegenes_counts_total$snp_eqtlgene),enriched_eqtls)
