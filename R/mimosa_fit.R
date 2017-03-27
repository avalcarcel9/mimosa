#' @title MIMoSA Fit on Training Data
#'
#' @description This function trains the MIMoSA model from a data.frame produced by an element from the output of the function mimosa_data.
#' @param training_dataframe data.frame(s) produced by the mimosa_data function
#' @param formula formula to be fit by glm model
#' @importFrom stats glm
#' @return Returns a glm object containing the trained MIMoSA coefficients.
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

