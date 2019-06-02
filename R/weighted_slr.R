#' @title Weighted simple linear regression
#'
#' @description Returns parameters from weighted regression of y on x (primary) and x on y (secondary, "reverse")
#' @param x independent modality matrix of nhood values (primary)
#' @param y dependent modality matrix of nhood values (primary)
#' @param wts vector of kernel weights
#' @param mInds indices for which to compute IMCo
#' @param refImg reference image for header info
#' @export
#' @return L list including intercept and slope from both weighted SLRs
# @examples \dontrun{
# 
#}
weighted_slr <- function(x, y, wts, mInds, refImg = NULL) {
  cs = function(x) {
    colSums(x, na.rm = TRUE)
  }
  # Missing data should be missing for both
  isna_x = is.na(x)
  isna_y = is.na(y)
  remove = isna_x | isna_y
  x[remove] = NA
  y[remove] = NA
  sumWtsMat = matrix(1, nrow=nrow(x), ncol=ncol(x))*wts
  sumWtsMat[remove] = NA
  rm(list = c("isna_x", "isna_y", "remove")); 
  gc()      
  # Compute weighted regression coefficient estimates using formulas
  sum_wts = cs(sumWtsMat)
  xw = x * wts
  yw = y * wts
  x2w = x^2 * wts
  y2w = y^2 * wts
  xwy = x * y * wts
  s_yw = cs(yw)
  s_xw = cs(xw)
  num = cs(xwy) - (s_yw * s_xw)/sum_wts
  denom = cs(x2w) - cs(xw)^2 / sum_wts
  denom_yx = cs(y2w) - cs(yw)^2 / sum_wts
  beta = num / denom
  beta_yx = num / denom_yx
  alpha = (s_yw - s_xw * beta) / sum_wts
  alpha_yx = (s_xw - s_yw * beta_yx) / sum_wts
  # Weighted R squared
  b1Mat = t(t(y)*beta)
  fitted = t(t(b1Mat) + alpha)
  resids = x - fitted
  m = t(t(fitted)/sum_wts)
  m = cs(wts*m)
  mss = cs(wts*t(t(fitted) - m)^2)
  rss = cs(wts*resids^2)
  tss = mss + rss
  r2 = mss/tss
  b1Mat_yx = t(t(x)*beta_yx)
  fitted_yx = t(t(b1Mat_yx) + alpha_yx)
  resids_yx = y - fitted_yx
  m_yx = t(t(fitted_yx)/sum_wts)
  m_yx = cs(wts*m_yx)
  mss_yx = cs(wts*t(t(fitted_yx) - m_yx)^2)
  rss_yx = cs(wts*resids_yx^2)
  tss_yx = mss_yx + rss_yx
  r2_yx = mss_yx/tss_yx
  # Return weighted regression coefficient estimates
  L = list(
    beta0 = alpha,
    beta1 = beta,
    r2 = r2,
    beta0_reverse = alpha_yx,
    beta1_reverse = beta_yx,
    r2_reverse = r2_yx)
  if (!is.null(refImg)) {
    L = lapply(L, make_ants_image,
               mask_indices = mInds,
               reference = refImg)
  }
  return(L)
}
