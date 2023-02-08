#install.packages("LDlinkR")

library('LDlinkR')

cells <- c("B", "CD4T", "CD8T", "NK", "DC", "monocyte")

output_file <- c()

# To run this scipt file with LD structure of given population is needed. After all the significant SNPs aremacthed with other SNPs in hight LD. ALternatively, we could also use LD package, which make the mapping of LD as well.  
outPath <- '/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/output/filtered_results/'

for (cell in cells){
  name <- paste(outPath, "UT_", cell, '/coeqtls_fullresults.sig.tsv.gz', sep="")
  tab <- read.table(name, sep='\t', header=T )
  # SNP <- as.data.frame(tab[,'SNP'])
  # colnames(SNP) <- "SNP"
  SNP <- (tab[,'SNP'])
  output_file <- c(SNP, output_file)
}

length(output_file)

output_file <- unique(output_file)
length(output_file)

LD_Score <- LDexpress(output_file[1], 
                                      pop = "CEU", 
                                      tissue = "ALL", 
                                      r2d = "r2", 
                                      r2d_threshold = 0.8, 
                                      p_threshold = 0.1, 
                                      win_size = 500000, 
                                      token = "d1bfc9a7a30b", 
                                      file = FALSE
)

for (i in output_file[63:length(output_file)] ){
print(i)
LD_Score_ind <- LDexpress(i, 
          pop = "CEU", 
          tissue = "ALL", 
          r2d = "r2", 
          r2d_threshold = 0.8, 
          p_threshold = 0.1, 
          win_size = 500000, 
          token = "d1bfc9a7a30b", 
          file = FALSE
)
LD_Score <- rbind(LD_Score_ind,LD_Score)
}

LD_Score_subset <- subset(LD_Score, select = c('Query',"RS_ID","R2"  ))
dim(LD_Score_subset)
LD_Score_subset <- LD_Score_subset[!duplicated(LD_Score_subset$RS_ID),]
dim(LD_Score_subset)


write.table(LD_Score,'/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/GRN_downstream_analysis/sign_LD_SNPs_18_12.txt', quote = F, col.names = F, row.names = F, sep='\t')

write.table(LD_Score_subset,'/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/GRN_downstream_analysis/sign_LD_SNPs_subset_18_12.txt', quote = F, col.names = F, row.names = F, sep='\t')

######### LDtrait



LD_Score <- LDtrait(output_file[2], 
                      pop = "CEU", 
                      r2d = "r2", 
                      r2d_threshold = 0.8, 
                      token = "d1bfc9a7a30b", 
                      file = FALSE
)
LD_Score

for (i in output_file){
  print(i)
      LD_Score_ind <- LDtrait(i, 
                              pop = "CEU", 
                              r2d = "r2", 
                              r2d_threshold = 0.8, 
                              token = "d1bfc9a7a30b", 
                              file = FALSE
      )
      LD_Score <- rbind(LD_Score_ind,LD_Score)
}

LD_Score_ind <- LDtrait(output_file[1:49], 
                        pop = "CEU", 
                        r2d = "r2", 
                        r2d_threshold = 0.8, 
                        token = "d1bfc9a7a30b", 
                        file = FALSE
)

LD_Score_ind2 <- LDtrait(output_file[50:72], 
                        pop = "CEU", 
                        r2d = "r2", 
                        r2d_threshold = 0.8, 
                        token = "d1bfc9a7a30b", 
                        file = FALSE
)
LD_Score <- rbind(LD_Score_ind,LD_Score_ind2)

# for (i in output_file){
#   print(i)
#   if(i %in% LD_Score$Query){
#     print('SNP is analyzed')
#   } else {
#     
#   tryCatch({
#     LD_Score_ind <- LDtrait(i, 
#                             pop = "CEU", 
#                             r2d = "r2", 
#                             r2d_threshold = 0.8, 
#                             token = "d1bfc9a7a30b", 
#                             file = FALSE
#     )
#     LD_Score <- rbind(LD_Score_ind,LD_Score)
#   }, error = function(e){
#     output_file <- output_file[-i]
#   print(length(output_file))
#   })
#   }
# }

# LD_Score_subset <- subset(LD_Score, select = c('Query',"RS_ID","R2"  ))
# dim(LD_Score_subset)
# LD_Score_subset <- LD_Score_subset[!duplicated(LD_Score_subset$RS_ID),]
# dim(LD_Score_subset)


write.table(LD_Score,'/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/GRN_downstream_analysis/sign_LD_SNPs_23_01.txt', quote = F, col.names = F, row.names = F, sep='\t')


#write.table(output_file,'/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/GRN_downstream_analysis/sign_SNPs_17_12.txt', quote = F, col.names = F, row.names = F, sep='\t')



expand_ld_table <- function(ld){
  # double the ld table, so we can easily select just from the left or right
  ld_copy <- ld[, c('CHR_B', 'BP_B', 'SNP_B', 'CHR_A', 'BP_A', 'SNP_A', 'R2')]
  colnames(ld_copy) <- c('CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B', 'SNP_B', 'R2')
  ld <- rbind(ld, ld_copy)
  # add each SNP in max LD with itself by copying the unique snps on the left and right
  ld_left <- ld[, c('CHR_A', 'BP_A', 'SNP_A')]
  ld_right <- ld[, c('CHR_B', 'BP_B', 'SNP_B')]
  colnames(ld_right) <- c('CHR_A', 'BP_A', 'SNP_A')
  ld_left_right <- rbind(ld_left, ld_right)
  ld_left_right <- unique(ld_left_right)
  ld_self <- cbind(ld_left_right, ld_left_right)
  colnames(ld_self) <- c('CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B', 'SNP_B')
  # ld with itself is off course 1
  ld_self$R2 <- 1
  # add to existing ld table
  ld <- rbind(ld, ld_self)
  return(ld)
}
# location of the LD file
ld_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LD_DB/genotypes_eur/EUR.chrAll.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.positions_plus_RSID.plink1.ldwindow10000.r2_075.ld'
# make the ld table a bit easier to work with
#ld_loc <- read.table(ld_loc, sep='\t')
ld <- expand_ld_table(ld_loc)
# confine the ld table to eQTL snps on the left <- subset here to what SNPs you need the LD of, with the other SNPs
ld <- ld[ld$SNP_A %in% eqtls$V1, ]