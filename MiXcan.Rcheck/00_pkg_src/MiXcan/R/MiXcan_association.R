#' Association analysis with MiXcan predicted genetically regulated expression
#'
#' The input MiXcan predicted genetically regulated expression levels are at cell-type
#' level. Association analysis returns p-value at tissue level for inference purpose, and
#' at cell-type level for post-inference (ad hoc) explaination of the likely source of
#' association.
#' @param new_y: A N by 2 matrix indicating the predicted gene expression levels in two
#' cell types. If a non-specific model is used for prediction, the predicted values should be the same in two cell types.
#' @param new_cov: A N by P matrix for the covariates in association analysis.
#' @param new_outcome: The dependent variable of the association analysis. It can be continuous, categorical
#' and count variables, and modeled through generalized linear models (GLM).
#' @param family: Same as the family object used in GLM.
#'
#' @return A tibble showing the estimate, standard error, p value for cell 1 and 2
#' and the combined p-value using the Cauchy combination.
#' @export
#'
MiXcan_association <- function(new_y, new_cov,
                               new_outcome, family= gaussian){
  dat <- data.frame(new_y, new_cov,  y = new_outcome)
  ft <- glm(as.formula(paste0("y ~.")), data = dat, family = family)
  res <- broom::tidy(ft)
  res1 <- res %>% filter(term == "cell_1")
  res2 <- res %>% filter(term == "cell_2")
  if( cor(new_y[,1],new_y[,2]) > 0.99) {res2 <- res1}
  if( cor(new_y[,1],new_y[,2]) < -0.99) {
    res2 <- res1 %>% mutate(estimate = -estimate)
  }
  p_combined <- ACAT::ACAT(res1$p.value, res2$p.value)
  result <- res1 %>%
    dplyr::select(cell1_est = estimate,
           cell1_se = std.error,
           cell1_p = p.value) %>%
    bind_cols(res2 %>%
                dplyr::select(cell2_est = estimate,
                       cell2_se = std.error,
                       cell2_p = p.value)) %>%
    mutate(p_combined = p_combined) %>%
    as.data.frame()
  return(result)
}
