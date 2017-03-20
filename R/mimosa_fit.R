#' @title Fits appropriate GLM model based on what image inputs are given
#'
#' @description Generates the MIMoSA model to create lesion probability maps, hard segementations in conjunction with mimosa_train_data
#' @param formula formula to be used on model
#' @param training_dataframe dataframe exported from mimosa_train_dataframe. list if training for multiple subjects
#' @param optimal_threshold a vector of potential thresholds to be used in the optimal thresholding algorithm
#' @param voxel_selection list of top voxels for subjects in training obtain from mimosa_train_dataframe or using top_voxel function in oasis package
#' @param gold_standard is the gold standard nifti object in list form
#' @export
#' @import fslr 
#' @import methods
#' @import nuerobase
#' @import oro.nifti
#' @import parallel
#' @import stats
#' @import oasis
#' @importFrom stats cov.wt qnorm
#' @return GLM objects fit in the MIMoSA procedure and optimal threshold
#' @examples \dontrun{
#' 
#'}

mimosa_fit <- function(training_dataframe, formula) {
  
  mimosa_model = glm(formula = formula, data = training_dataframe, 
                      family = binomial)
  mimosa_model$y = c()
  mimosa_model$model = c()
  mimosa_model$residuals = c()
  mimosa_model$fitted.values = c()
  mimosa_model$effects = c()
  mimosa_model$qr$qr = c()
  mimosa_model$linear.predictors = c()
  mimosa_model$weights = c()
  mimosa_model$prior.weights = c()
  mimosa_model$data = c()
  mimosa_model$family$variance = c()
  mimosa_model$family$dev.resids = c()
  mimosa_model$family$aic = c()
  mimosa_model$family$validmu = c()
  mimosa_model$family$simulate = c()
  attr(mimosa_model$terms, ".Environment") = c()
  attr(mimosa_model$formula, ".Environment") = c()
  
  return(mimosa_model)
}

