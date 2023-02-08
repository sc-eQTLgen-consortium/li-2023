# ------------------------------------------------------------------------------
# Generate an upset plot of overlap between cell types
# Input: significant co-eQTL results per cell type
# Output: upset plot
# ------------------------------------------------------------------------------

library(data.table)
library(UpSetR)
  
coeqtls_mono<-fread("coeqtl_mapping/output/filtered_results/UT_monocyte/coeqtls_fullresults_fixed.sig.tsv.gz")
coeqtls_cd4t<-fread("coeqtl_mapping/output/filtered_results/UT_CD4T/coeqtls_fullresults_fixed.sig.tsv.gz")
coeqtls_cd8t<-fread("coeqtl_mapping/output/filtered_results/UT_CD8T/coeqtls_fullresults_fixed.sig.tsv.gz")
coeqtls_nk<-fread("coeqtl_mapping/output/filtered_results/UT_NK/coeqtls_fullresults_fixed.sig.tsv.gz")
coeqtls_dc<-fread("coeqtl_mapping/output/filtered_results/UT_DC/coeqtls_fullresults_fixed.sig.tsv.gz")
coeqtls_b<-fread("coeqtl_mapping/output/filtered_results/UT_B/coeqtls_fullresults_fixed.sig.tsv.gz")

pdf(paste0(outdir, "grn_plot_snp_gene_gene.pdf"))

upset(fromList(list(Monocyte = coeqtls_mono$snp_genepair,
                    `CD4+ T` = coeqtls_cd4t$snp_genepair,
                    `CD8+ T` = coeqtls_cd8t$snp_genepair,
                    NK = coeqtls_nk$snp_genepair,
                    DC = coeqtls_dc$snp_genepair,
                    B = coeqtls_b$snp_genepair)), 
      set_size.show = T,set_size.scale_max = 600, 
      mainbar.y.label = "SNP-Gene-Gene", 
      nintersects = 40, nsets = 10,
      text.scale = 1.5)

dev.off()

#Identify all elements that are in at least four of the six cell types
all_coeqtls<-c(unique(coeqtls_mono$snp_genepair),
               unique(coeqtls_cd4t$snp_genepair),
               unique(coeqtls_cd8t$snp_genepair),
               unique(coeqtls_nk$snp_genepair),
               unique(coeqtls_dc$snp_genepair),
               unique(coeqtls_b$snp_genepair))

occurrence<-data.frame(table(all_coeqtls))

#Show all coeQTLs part of at least three different cell types:
most_occ<-occurrence[occurrence$Freq > 2,]
#How many of the frequent coeQTls are associated with the RPS26 locus:
mean(startsWith(as.character(most_occ$all_coeqtls),"rs1131017"))

