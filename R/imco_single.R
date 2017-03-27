#' @title Single Inter-modal Coupling model
#'
#' @description Data for coupling estimation on a single subject and single voxel
#' @param files Vector of full paths to images
#' @param vxl Vector of single voxel indicies to get data from
#' @param brainMask Full path to brain mask
#' @param ref reference modality or dependent modality when type="regression"
#' @param neighborhoodSize Full width in voxels for FWHM
#' @import ANTsR
#' @import extrantsr
#' @import rlist
#' @return Data frame
#' @examples \dontrun{
#'
#'}
imco_single <- function(files, vxl, brainMask, ref=1, neighborhoodSize=3){
  if(!(as.integer(neighborhoodSize)==neighborhoodSize)){
    stop('neighborhoodSize must be an odd integer')
  }
  if(!(is.numeric(ref) & ref<=length(files))){
    stop('check reference modality specification')
  }
  nf = length(files)
  fileList = check_ants(files)
  for(i in 2:length(files)){
    if(!all(dim(fileList[[i-1]])==dim(fileList[[i]]))){
      stop('Image dimensions do not match')
    }
    if(!all(antsGetSpacing(fileList[[i-1]])==antsGetSpacing(fileList[[i]]))){
      stop('Voxel dimensions do not match')
    }
  }
  # Dimension of each voxel in mm
  vDims = antsGetSpacing(fileList[[1]])
  # Image dimension
  imgDims = dim(fileList[[1]])
  # Read in brain mask
  bMask = check_ants(brainMask)
  if(!all(dim(bMask)==dim(fileList[[1]]))){
    stop('Image dimensions do not match the bran mask')
  }
  if(!all(antsGetSpacing(bMask)==antsGetSpacing(fileList[[1]]))){
    stop('Voxel dimensions do not match the brain mask')
  }
  # Assign anything outside brain mask to NA
  fileList = lapply(fileList, function(x){newx = x; newx[bMask==0] = NA; return(newx)})
  ######################################
  # FWHM => sigma
  # Note: We specify FWHM in terms of number of
  #       voxels along the axis of smallest size (mm).
  #       Thus, the neighborhood is specified in terms
  #       of a univariate Gaussian but the weights are
  #       computed from a trivariate Gaussian consisting
  #       of 3 independent univariate Gaussians.
  ######################################
  minDim = min(vDims) #Should I use max or min to compute sigma?
  if(neighborhoodSize%%2!=1){
    stop('neighborhoodSize must be odd integer')
  }
  width = neighborhoodSize*minDim
  # Need sigma for specifying Guassian kernel
  sigma = width/(2*sqrt(2*log(2)))
  # We will choose a radius based on quantile (mm) of univariate normal
  lower = qnorm(.005, sd=sigma)
  # Scale by voxel dimensions
  radius = floor(abs(lower)/vDims)
  # Neighborhood data from each modality
  if(!is.null(subMask)){
    sMask = check_ants(subMask)
    nhoods = lapply(fileList, function(x) getNeighborhoodInMask(image=x, mask=sMask, radius=radius, spatial.info=TRUE, boundary.condition='image'))
  } else{
    nhoods = lapply(fileList, function(x) getNeighborhoodInMask(image=x, mask=bMask, radius=radius, spatial.info=TRUE))
  }
  # Will use to map back to niftis
  inds = nhoods[[1]]$indices
  inds1 = inds + 1
  # Will use to compute distances from center voxel
  offs = nhoods[[1]]$offsets
  dists = t(apply(offs, 1, function(x, y) x*y, y=vDims))
  sqNorm = apply(dists*dists, 1, function(x) sum(x))
  # Constant doesn't matter because we would rescale by
  # max weight so that center voxel would be be 1
  nWts = exp(-1*sqNorm/(2*sigma^2))
  imgVals = lapply(nhoods, function(x) x$values)
  vInd = apply(inds1, 1, function(x) all(vxl==x))
  getData = list.cbind(lapply(imgVals, function(x) x[,vInd]))
  getData = cbind(getData[,ref], getData[,-ref])
  return(list('X'=getData, 'w'=nWts))
}
