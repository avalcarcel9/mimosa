#' @title Inter-modal Coupling Model
#'
#' @description Implements full volumetric IMCo coupling estimation on a single subject
#' @param files Vector of full paths to images or a list of local nifti objects
#' @param brainMask Full path to brain mask or local nifti object
#' @param subMask NULL or local nifti object or the full path to a mask that is a subset of brain mask where coupling should be computed
#' @param type "regression" or "pca"
#' @param ref reference modality when type="pca" or dependent modality when type="regression"
#' @param neighborhoodSize Full width in voxels for FWHM
#' @param reverse If TRUE, calculates both regressions if type="regression" and length(files)==2, otherwise ignored
#' @param verbose TRUE for updates on computation, else FALSE
#' @param retimg If TRUE, return list of estimated coupling maps as nifti objects
#' @param outDir Full path to directory where maps should be written
#' @importFrom ANTsRCore antsGetSpacing getNeighborhoodInMask
#' @importFrom extrantsr check_ants
#' @importFrom stats qnorm
#' @return Estimated IMCo coupling maps, either written to files and/or returned as nifti objects
#' @examples \dontrun{
#'
#'}
imco <- function(files, brainMask, subMask=NULL, type="pca", ref=1, neighborhoodSize=3, reverse=TRUE, verbose=TRUE, retimg=FALSE, outDir=NULL){
    if(!(type=="pca" | type=="regression")){
        stop('type must be either pca or regression')
    }
    if(!(is.numeric(ref) & ref<=length(files))){
        stop('check reference modality specification')
    }
    if(!(as.integer(neighborhoodSize)==neighborhoodSize)){
        stop('neighborhoodSize must be an odd integer')
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
        stop('Image dimensions do not match the brain mask')
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
    if(verbose){
    	cat("# Extracting neighborhood data \n")
    }
    # Neighborhood data from each modality
    if(!is.null(subMask)){
        sMask = check_ants(subMask)
        mask_indices = which(as.array(sMask) > 0)
        nhoods = lapply(fileList, function(x) getNeighborhoodInMask(image=x, mask=sMask, radius=radius, spatial.info=TRUE, boundary.condition='image'))
    } else{
        mask_indices = which(as.array(bMask) > 0)
        nhoods = lapply(fileList, function(x) getNeighborhoodInMask(image=x, mask=bMask, radius=radius, spatial.info=TRUE))
    }
    # Will use to map back to niftis
    inds = nhoods[[1]]$indices
    offs = nhoods[[1]]$offsets
    # ANTsR version 0.3.2 used 0 first index instead of 1
    # This was changed here: https://github.com/stnava/ANTsR/commit/706148aa994d414c9efd76e9292ae99351e2c4be
    # Thus, have commented out the following line and require latest release of ANTsR version 0.3.3
    # inds = inds + 1
    # Will use to compute distances from center voxel
    nWts = get_weights(offs, vDims, sigma=sigma)
    if(type=="regression"){
    	regObj = imco_reg(files=files, nhoods=nhoods, nWts=nWts, mask_indices=mask_indices, ref=ref, reverse=reverse, verbose=verbose, retimg=retimg, outDir=outDir)
    	return(regObj)
    } else{
    	pcaObj = imco_pca(files=files, nhoods=nhoods, nWts=nWts, mask_indices=mask_indices, ref=ref, verbose=verbose, retimg=retimg, outDir=outDir)
    	return(pcaObj)
    }
}


