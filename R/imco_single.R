#' @title Single Inter-modal Coupling model
#'
#' @description Data for coupling estimation on a single subject and single voxel
#' @param files Vector of full paths to images or a list of local nifti objects
#' @param vxl Vector of single voxel indicies to get data from
#' @param brainMask Full path to brain mask or local nifti object
#' @param subMask NULL or local nifti object or the full path to a mask that is a subset of brain mask where coupling should be computed
#' @param ref reference modality when type="pca" or dependent modality when type="regression"
#' @param neighborhoodSize Full width in voxels for FWHM
#' @export
#' @importFrom extrantsr check_ants
#' @importFrom rlist list.cbind
#' @return Data frame
#' @examples \dontrun{
#' 
#'}
imco_single <- function(files, vxl, brainMask, subMask=NULL, ref=1, neighborhoodSize=3){
    if(!(as.integer(neighborhoodSize)==neighborhoodSize)){
        stop('neighborhoodSize must be an odd integer')
    }
    if(!(is.numeric(ref) & ref<=length(files))){
        stop('check reference modality specification')
    }
    nf = length(files)
    fileList = extrantsr::check_ants(files)
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
    # Will use to map back to niftis - these are the center indices, one for each voxel in the gray matter mask
    inds = nhoods[[1]]$indices
    # ANTsR version 0.3.2 used 0 first index instead of 1
    # This was changed here: https://github.com/stnava/ANTsR/commit/706148aa994d414c9efd76e9292ae99351e2c4be
    # Thus, have commented out the following line and require latest release of ANTsR version 0.3.3
    # inds = inds + 1
    # Will use to compute distances from center voxel
    offs = nhoods[[1]]$offsets
    dists = t(apply(offs, 1, function(x, y) x*y, y=vDims))
    sqNorm = apply(dists*dists, 1, function(x) sum(x))
    # Constant doesn't matter because we would rescale by
    # max weight so that center voxel would be be 1
    nWts = exp(-1*sqNorm/(2*sigma^2))
    imgVals = lapply(nhoods, function(x) x$values)
    vInd = apply(inds, 1, function(x) all(vxl==x))
    getData = list.cbind(lapply(imgVals, function(x) x[,vInd]))
    getData = cbind(getData[,ref], getData[,-ref])
    return(list('X'=getData, 'w'=nWts))
}

