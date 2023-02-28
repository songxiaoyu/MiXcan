#' Predict the cell-type specific or non-specific expression levels of a gene with MiXcan model in new genetic data.
#'
#' @param weight Weight matrix produced by MiXcan_extract_weight()
#' @param new_x The genetic data used for prediction in a new data set.
#'
#' @return A N by 2 matrix indicating the predicted gene expression levels in two
#' cell types. If a non-specific model is used for prediction, the predicted values should be the same in two cell types.
#' @export
#'
MiXcan_prediction <- function(weight, new_x){
  yhat_MiXcan_cell_1 <- new_x %*% as.matrix(weight[,"weight_cell_1"])
  yhat_MiXcan_cell_2 <- new_x %*% as.matrix(weight[,"weight_cell_2"])
  yhat_MiXcan_prediction <- cbind(yhat_MiXcan_cell_1, yhat_MiXcan_cell_1)
  colnames(yhat_MiXcan_prediction) <- c("cell_1", "cell_2")
  return(yhat_MiXcan_prediction)

}
