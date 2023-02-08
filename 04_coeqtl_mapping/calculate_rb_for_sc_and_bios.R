# Title     : TODO
# Objective : TODO
# Created by: Shuang
# Created on: 1/19/2022

source("Rb.R")
library(glue)
library(data.table)

calculate_rb_bios_replication_summary <- function(biostype, filtertype){
    print(biostype)
    print(filtertype)
    resdf <- c()
    for ( celltype in c('CD4T', 'CD8T', 'monocyte', 'NK', 'B', 'DC') ){
        df <- read.csv(glue('./coeqtl_mapping/bios/{biostype}/{filtertype}/UT_{celltype}/replication_parameters.csv'))
        res <- calcu_cor_true(df$flipped_bios_beta, df$std.err_bios, df$MetaBeta, df$MetaSE, df$theta)
        res <- cbind(res, celltype)
        resdf <- rbind(res, resdf)
    }
    write.csv(resdf, glue('./coeqtl_mapping/bios/{biostype}/{filtertype}/replication_summary.csv'))
}

# BIOS replication
args = commandArgs(trailingOnly=TRUE)
calculate_rb_bios_replication_summary(args[1], args[2])


# coeQTLs
filtertype = 'filtered_results'
workdir = './coeqtl_mapping/'
resdf <- c()
for ( celltype_discovery in c('CD4T', 'CD8T', 'monocyte', 'NK', 'B', 'DC') ){
    for ( celltype_replication in c('CD4T', 'CD8T', 'monocyte', 'NK', 'B', 'DC') ){
        if ( celltype_discovery != celltype_replication ){
            df <- fread(glue('{workdir}/output/{filtertype}/rb_calculations/discovery_{celltype_discovery}_replication_{celltype_replication}.tsv.gz'))
            print(c(celltype_discovery, celltype_replication, nrow(df)))
            if ( nrow(df) < 5 ){
                resdf <- rbind(resdf, c(NA, NA, 0, celltype_discovery, celltype_replication))
            }else{
                res <- calcu_cor_true(df$MetaBeta, df$MetaSE, df$MetaBeta_replication, df$MetaSE_replication, df$theta)
                res <- cbind(res, celltype_discovery, celltype_replication)
                resdf <- rbind(res, resdf)
            }
        }
    }
}
write.csv(resdf, glue('{workdir}/output/{filtertype}/rb_calculations/summary.csv'))

# coeQTLs monocyte sub celltypes
filtertype = 'filtered_results'
workdir = './coeqtl_mapping/'
resdf <- c()
for ( celltype_discovery in c('monocyte', 'cMono', 'ncMono') ){
    for ( celltype_replication in c('monocyte', 'cMono', 'ncMono') ){
        if ( celltype_discovery != celltype_replication ){
            df <- fread(glue('{workdir}/output/{filtertype}/rb_calculations/monocyte_subcelltypes/discovery_{celltype_discovery}_replication_{celltype_replication}.tsv.gz'))
            print(c(celltype_discovery, celltype_replication, nrow(df)))
            if ( nrow(df) < 5 ){
                resdf <- rbind(resdf, c(NA, NA, 0, celltype_discovery, celltype_replication))
            }else{
                res <- calcu_cor_true(df$MetaBeta, df$MetaSE, df$MetaBeta_replication, df$MetaSE_replication, df$theta)
                res <- cbind(res, celltype_discovery, celltype_replication)
                resdf <- rbind(res, resdf)
            }
        }
    }
}
write.csv(resdf, glue('{workdir}/output/{filtertype}/rb_calculations/monocyte_subcelltypes/summary.csv'))


# eQTLs
workdir = './coeqtl_mapping/'
resdf <- c()
for ( celltype_discovery in c('CD4T', 'CD8T', 'monocyte', 'NK', 'B', 'DC') ){
    for ( celltype_replication in c('CD4T', 'CD8T', 'monocyte', 'NK', 'B', 'DC') ){
        if ( celltype_discovery != celltype_replication ){
            df <- fread(glue('{workdir}/input/snp_selection/rb_calculations/discovery_{celltype_discovery}_replication_{celltype_replication}.tsv.gz'))
            print(c(celltype_discovery, celltype_replication, nrow(df)))
            if ( nrow(df) < 5 ){
                resdf <- rbind(resdf, c(NA, NA, 0, celltype_discovery, celltype_replication))
            }else{
                res <- calcu_cor_true(df$metabeta, df$SE, df$metabeta_replication, df$SE_replication, df$theta)
                res <- cbind(res, celltype_discovery, celltype_replication)
                resdf <- rbind(res, resdf)
            }
        }
    }
}
write.csv(resdf, glue('{workdir}/input/snp_selection/rb_calculations/summary.csv'))