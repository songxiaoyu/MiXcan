#' Extract the gene expression prediction summary
#'
#' Extract the gene expression prediction summary (e.g. No of SNPs, R^2, correlation) function.
#' The output can be directly applied to GWAS data for cell-type specific TWAS.
#'
#' @param model A direct output from MiXcan() function, which includes the
#' prediction weights as well as other summaries of the prediction models.
#' @param y: The pre-cleaned expression level data for a single gene in N samples.
#' @param x: A N by P matrix for all the genetic predictors used to predict the genetically regulated expression  of the gene.
#' @param pi: An estimation of cell-type faction of the cell type of interest.
#'
#' @return A data frame with weight for cell 1 and 2, including the potential meta data for the SNP
#' @export
#'
MiXcan_extract_summary <- function(model, x, y, pi) {

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


#' Extract the gene expression prediction summary
#'
#' Extract the gene expression prediction summary (e.g. No of SNPs, R^2, correlation) function.
#' The output can be directly applied to GWAS data for cell-type specific TWAS.
#'
#' @param model A direct output from MiXcan() function, which includes the
#' prediction weights as well as other summaries of the prediction models.
#' @param y: The pre-cleaned expression level data for a single gene in N samples.
#' @param x: A N by P matrix for all the genetic predictors used to predict the genetically regulated expression  of the gene.
#' @param cov: A N by Q matrix for the covariates adjusted in the model (e.g. age, population stratification).
#' @param pi: An estimation of cell-type fraction for the cell type of interest.
#' 

#' @return A data frame with weight for cell 1 and 2, including the potential meta data for the SNP
#' @export
#'
MiXcan_refit_summary <- function(model, x, y, pi, cov=NULL) {
  summary <- MiXcan_extract_summary(x=x, y=y, pi=pi, model=model)
  weights=MiXcan_extract_weight(model, keepZeroWeight=T)
  w2 <- MiXcan_refit_weight(model = model, y=y,x=x, cov = NULL, pi= pi)
  idx=match(w2$xNameMatrix, weights$xNameMatrix)
  weights$weight_cell_1[idx]=w2$weight_cell_1 
  weights$weight_cell_2[idx]=w2$weight_cell_2 
  
  yhat=x %*% weights$weight_cell_1 *pi + x %*% weights$weight_cell_2 *(1-pi)
  
  r=cor.test(y,yhat, use="complete.obs")
  summary2= summary %>% data.frame() %>% 
    mutate(in_sample_r2_refit=r$estimate^2) %>%
    mutate(in_sample_cor_pvalue_refit = r$p.value)

  return(summary2)
}






