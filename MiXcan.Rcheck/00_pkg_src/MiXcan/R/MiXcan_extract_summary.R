#' Extract the gene expression prediction weights matrix
#' 
#' Extract the gene expression prediction weights matrix (elastic net regression coefficients) from MiXcan() function.
#' The output can be directly applied to GWAS data for cell-type specific TWAS.
#'
#' @param MiXcan_model A direct output from MiXcan() function, which includes the
#' prediction weights as well as other summaries of the prediction models.
#'
#' @return A data frame with weight for cell 1 and 2, including the potential meta data for the SNP
#' @export
#'
MiXcan_extract_summary <- function(x, y, pi, foldid, model) {
  
  weights=MiXcan_extract_weight(model, keepZeroWeight=T)
  w2=MiXcan_extract_weight(model, keepZeroWeight=F)
  
  yhat=x %*% weights$weight_cell_1 *pi + x %*% weights$weight_cell_2 *(1-pi)
  
  r=cor.test(y,yhat, use="complete.obs")
  
  summary=cbind(model$yName, n_snp_input=ncol(x),
                n_snp_model=nrow(w2),
                model_type=model$type,
                in_sample_r2=r$estimate^2,
                in_sample_cor_pvalue=r$p.value)
  
  return(summary)
}



