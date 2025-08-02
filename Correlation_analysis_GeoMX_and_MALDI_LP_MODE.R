# The general code for spatial transcriptomics data and correlation with MALDI MSI LP mode FWHM:

#corrFunc <- function(var1, var2, data1) {
#  if (!is.numeric(data1[[var1]]) || !is.numeric(data1[[var2]])) {
#    warning("Both variables must be numeric for correlation calculation.")
#    return(NULL)
#  }
  
#  result <- cor.test(data1[[var1]], data1[[var2]], na.rm = TRUE)
#  data.frame(var1, var2, result[c("estimate", "p.value", "statistic", "method")], 
#             stringsAsFactors = FALSE)
#}


# corrs_trans_df_AB_for_correlation_norm_Quant_10 = do.call(rbind, mapply(corrFunc, vars[,1], vars[,2], 
#                                                                        MoreArgs=list(data=trans_df_AB_for_correlation_norm_Quant_10), SIMPLIFY=FALSE))

# corrs_trans_df_AB_for_correlation_norm_Quant_18 = do.call(rbind, mapply(corrFunc, vars[,1], vars[,2], 
#                                                                        MoreArgs=list(data=trans_df_AB_for_correlation_norm_Quant_18), SIMPLIFY=FALSE))

