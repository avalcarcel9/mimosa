#' MIMoSA Model Formula
#'
#' @param T2 Logical. Should the terms including \code{T2} be included?
#' @param PD Logical. Should the terms including \code{PD} be included?
#'
#' @return A \code{formula} object 
#' @export
#'
#' @importFrom stats as.formula
#' @examples
#' mimosa_formula()
#' mimosa_formula(T2 = FALSE)
#' mimosa_formula(PD = FALSE)
#' mimosa_formula(T2 = FALSE, PD = FALSE)
mimosa_formula = function(T2 = TRUE, PD = TRUE) {
  formula = gold_standard ~ FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + 
    PD_10 * PD + PD_20 * PD + 
    T2_10 * T2 + T2_20 * T2 + 
    T1_10 * T1 + T1_20 * T1 +
    FLAIRonT1_intercepts + FLAIRonT2_intercepts + FLAIRonPD_intercepts +
    T1onT2_intercepts + T1onPD_intercepts + T2onPD_intercepts +
    T1onFLAIR_intercepts + T2onFLAIR_intercepts + PDonFLAIR_intercepts + 
    T2onT1_intercepts + PDonT1_intercepts + PDonT2_intercepts +
    FLAIRonT1_slopes + FLAIRonT2_slopes + FLAIRonPD_slopes +
    T1onT2_slopes + T1onPD_slopes + T2onPD_slopes +
    T1onFLAIR_slopes + T2onFLAIR_slopes + PDonFLAIR_slopes +
    T2onT1_slopes + PDonT1_slopes + PDonT2_slopes
  formula = as.character(formula)
  yvar = formula[2]
  formula = formula[3]
  run_terms = strsplit(formula, split = "+", fixed = TRUE)[[1]]
  run_terms = trimws(run_terms)
  if (!T2) {
    run_terms = run_terms[ !grepl("T2", run_terms)]
  }
  if (!PD) {
    run_terms = run_terms[ !grepl("PD", run_terms)]
  }  
  run_terms = paste(run_terms, collapse = " + ")
  formula = paste0(yvar, " ~ ", run_terms)
  formula = as.formula(formula)
  return(formula)
}
