# ------------------------------------------------------------------------------
# Fit logarithmic curve describing mean correlation between individuals dependent
# on the number of cells per individual and cell type, fit down separately
# for each cell type
# Input: pairwise comparison of all individuals (Pearson correlation) per cell type
#        and for different numbers of cells 
#        (subsampling and calculation done in correlation_subsampling.py)
# Output: logarithmic fit per cell type and curve visualizing the fit
# ------------------------------------------------------------------------------

library(ggplot2)
library(dplyr)

theme_set(theme_bw())

suffix<-"v3"
#suffix<-"v2"

color_coding <- list()
color_coding[["CD4+ T"]] <- "#2E9D33"
color_coding[["CD8+ T"]] <- "#126725"
color_coding[["Monocyte"]] <- "#EDBA1B"
color_coding[["NK"]] <- "#E64B50"
color_coding[["B"]] <- "#009DDB"
color_coding[["DC"]] <- "#965EC8"

#Full cell type names as reported in the paper
cell_types_corrected<-setNames(c("CD4+ T","CD8+ T","Monocyte","NK","DC","B"),
                               c("CD4T","CD8T","monocyte","NK","DC","B"))

#Load results
res<-read.csv(paste0("co-expression_indivs_subsampled/",
                     "correlation_individuals_subsampled_1M_",suffix,".csv"),
              stringsAsFactors = FALSE)
res$X<-NULL

res$celltype<-cell_types_corrected[res$celltype]

res_summary<-res%>%
  group_by(celltype,cell_num)%>%
  summarise(mean_corr=mean(corr),
            quantile_25=quantile(corr,0.25),
            quantile_75=quantile(corr,0.75))%>%
  as.data.frame()

#Filter out B cells and DCs because no line can be drawn for them
res_summary<-res_summary[! res_summary$celltype %in% c("B","DC"),]

#Fit one log function for each cell type
log_parameters<-NULL
for(cell_type in unique(res_summary$celltype)){
  
  #Fit the linear model
  res_ct<-res[res$celltype == cell_type,]
  model_lm<-lm(corr~log(cell_num),data=res_ct)
  
  #Save model summary
  summary_model<-summary(model_lm)
  print(summary_model)
  
  log_parameters<-rbind(log_parameters,
                        data.frame(cell_type,
                                   intercept=summary_model$coefficients[1,1],
                                   log_beta=summary_model$coefficients[2,1],
                                   adj_r_squared=summary_model$adj.r.squared,
                                   stringsAsFactors = FALSE))
  res_summary<-rbind(res_summary,
                     data.frame(celltype=cell_type,
                                cell_num=seq(max(res_ct$cell_num)+25,1500,by=25),
                                mean_corr=NA,
                                quantile_25=NA,
                                quantile_75=NA))
}



res_summary$fitted_corr<-sapply(1:nrow(res_summary),function(i)
  log_parameters$intercept[log_parameters$cell_type == res_summary$celltype[i]] +
  log_parameters$log_beta[log_parameters$cell_type == res_summary$celltype[i]] *
  log(res_summary$cell_num[i]))

res_summary_melt<-reshape2::melt(res_summary[,c("celltype","cell_num","mean_corr","fitted_corr")],
                                 id.vars=c("celltype","cell_num"))

g<-ggplot()+
  geom_line(data=res_summary_melt,aes(x=cell_num,y=value,color=celltype,
                                      linetype=variable))+
  geom_point(data=res_summary,aes(x=cell_num,y=mean_corr,color=celltype))+
  # annotate("text",x=850,y=0.2,hjust=0,
  #          label=paste0("y ~ -0.56 + 0.21 * log(x), R^2 = 0.98 (CD4+ T)\n",
  #                       "y ~ -0.48 + 0.20 * log(x), R^2 = 0.86 (CD8+ T)\n",
  #                       "y ~ -0.53 + 0.20 * log(x), R^2 = 0.94 (Monocyte)\n",
  #                       "y ~ -0.41 + 0.15 * log(x), R^2 = 0.93 (NK)\n"))+
  scale_color_manual("Cell type",values=unlist(color_coding))+
  scale_linetype_discrete("",labels=c("Observed","Predicted"))+
  xlab("Subsampled number of cells per individual")+
  ylab("Correlation between individuals")
print(g)

ggsave(g,file=paste0("co-expression_indivs_subsampled/plots/subsampling_1M_",
                     suffix,"_fitted_lines.pdf"),
       width=10,height=4)
