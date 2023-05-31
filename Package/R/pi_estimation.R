#' A function to get robust cell-type composition estimation based on the prior estimates and transcriptomic data.
#'
#' @param expression_matrix: G by N gene expression matrix for estimating cell type composition of the tissue.
#' Each row is a gene and each column is tissue sample. Gene name and sample ID can be included as row and column names,
#' but not in the expression matrix.
#' @param n_iteration: Number of bootstrap samples to estimate the pi.
#' @param prior: Prior estimation of cell-type composition.
#'
#' @return Updated cell-type composition estimates. A ribble of 2 columns: The first columns shows the sample id, the second column shows the robust cell-type composition estimation
#' @export
#'
#'

pi_estimation <- function(expression_matrix,
                          n_iteration=5, prior){

  TSNetB_prop=foreach (i = 1:n_iteration, .combine = "rbind") %dopar%{
    sample_index <- sample(sample(1:dim(expression_matrix)[2], round(dim(expression_matrix)[2] * 0.8), replace = FALSE))
    TSNetB=try(deNet_purity(t(expression_matrix[,sample_index]),purity=prior[sample_index]))

    if(!inherits(TSNetB, "try-error")) {
      TSNetB_prop_once <- tibble(sample = colnames(expression_matrix[,sample_index]),
                                 prop = TSNetB[[1]],
                                 rep = i)
      TSNetB_prop_once
    }
  }

  TSNetB_prop_trim_mean <- TSNetB_prop %>%
    group_by(sample) %>%
    summarise(mean_trim_0.05 = mean(prop, trim = 0.05)) %>%
    ungroup
  return(TSNetB_prop_trim_mean)
}


