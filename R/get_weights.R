#' @title Weights for IMCo model
#'
#' @description Gets weights for IMCo model
#' @param offsets matrix of offsets from center voxel that define neighborhood
#' @param voxelDims voxel dimensions in mm
#' @param sigma standard deviation of Gaussian kernel
#' @export
#' @return vector of weights
#' @examples \dontrun{
#' 
#'}

get_weights <- function(offsets, voxelDims, sigma){
    dists = t(apply(offsets, 1, function(x, y) x*y, y=voxelDims))
    sqNorm = apply(dists*dists, 1, function(x) sum(x))
    # Constant doesn't matter because we would rescale by
    # max weight so that center voxel would be be 1
    return(exp(-1*sqNorm/(2*sigma^2)))
}