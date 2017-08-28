#' @title Weights for IMCo model
#'
#' @description Gets weights for IMCo model
#' @param offsets matrix of offsets from center voxel that define neighborhood
#' @param voxelDims voxel dimensions in mm
#' @param sigma standard deviation of Gaussian kernel
#' @export
#' @return vector of weights
#' @examples
#' offsets = c(-1, 0, 1)
#' offsets = expand.grid(offsets, offsets, offsets)
#' sigma = 3
#' voxelDims = c(0.25, 0.8, 3)
#' wts = get_weights(offsets = offsets, voxelDims = voxelDims
#' sigma = sigma)
#' stopifnot(length(wts)==27)
get_weights <- function(offsets, voxelDims, sigma){
  ndim = length(voxelDims)
  transpose = ncol(offsets) == ndim
  dont_transpose = nrow(offsets) == ndim
  if (!transpose && !dont_transpose) {
    stop("Matrix of offsets are not correct!")
  }
  if (transpose) {
    offsets = t(offsets)
  }
  offsets = offsets * voxelDims
  sqNorm = colSums(offsets^2)

  # dists = t(apply(offsets, 1, function(x, y) x*y, y=voxelDims))
  # sqNorm = apply(dists*dists, 1, function(x) sum(x))
  # Constant doesn't matter because we would rescale by
  # max weight so that center voxel would be be 1
  return(exp(-1*sqNorm/(2*sigma^2)))
}
