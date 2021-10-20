#' Association analysis with MiXcan predicted gene expression levels
#'
#' @param MiXcan_predicted_expr: A N by 2 matrix indicating the predicted gene expression levels in two
#' cell types. If a non-specific model is used for prediction, the predicted values should be the same in two cell types.
#' @param covariates: A N by P matrix for the covariates in association analysis.
#' @param outcome: The dependent variable of the association analysis. It can be continuous, categorical
#' and count variables, and modeled through generalized linear models (GLM).
#' @param family: Same as the family object used in GLM.
#'
#' @return A tibble showing the estimate, standard error, p value for cell 1 and 2 and the combined p value using the Cauchy combination. 
#' @export
#'
#' @examples
MiXcan_association <- function(MiXcan_predicted_expr, covariates, outcome, family){
  MiXcan_association_dataframe <- cbind(MiXcan_prediction_result, covariates,  y = outcome)
  MiXcan_association_dataframe_glm <- glm(as.formula(paste0("y ~.")), data = MiXcan_association_dataframe, family = family)
  MiXcan_association_dataframe_glm_result <- broom::tidy(MiXcan_association_dataframe_glm)
  MiXcan_association_dataframe_glm_result_cell_1 <- MiXcan_association_dataframe_glm_result %>%
    filter(term == "cell_1")
  MiXcan_association_dataframe_glm_result_cell_2 <- MiXcan_association_dataframe_glm_result %>%
    filter(term == "cell_2")
  if( cor(MiXcan_predicted_expr[,1],MiXcan_predicted_expr[,2]) > 0.99) {
    MiXcan_association_dataframe_glm_result_cell_2 <- MiXcan_association_dataframe_glm_result_cell_1
  }
  if( cor(MiXcan_predicted_expr[,1],MiXcan_predicted_expr[,2]) < -0.99) {
    MiXcan_association_dataframe_glm_result_cell_2 <- MiXcan_association_dataframe_glm_result_cell_1 %>%
      mutate(estimate = -estimate)
  }
  p_combined <- ACAT::ACAT(MiXcan_association_dataframe_glm_result_cell_1$p.value, MiXcan_association_dataframe_glm_result_cell_2$p.value)
  MiXcan_association_dataframe_result <- MiXcan_association_dataframe_glm_result_cell_1 %>%
    dplyr::select(cell_1_estimate = estimate,
           cell_1_std.error = std.error,
           cell_1_p = p.value) %>%
    bind_cols(MiXcan_association_dataframe_glm_result_cell_2 %>%
                dplyr::select(cell_2_estimate = estimate,
                       cell_2_std.error = std.error,
                       cell_2_p = p.value)) %>%
    mutate(p_combined = p_combined) %>%
    as.data.frame()
  return(MiXcan_association_dataframe_result)
}
