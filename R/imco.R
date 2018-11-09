#' @title Inter-modal Coupling Model
#'
#' @description Implements full volumetric IMCo coupling estimation on a single subject
#' @param files Vector of full paths to images or a list of local nifti objects
#' @param brainMask Full path to brain mask or local nifti object
#' @param subMask NULL or local nifti object or the full path to a mask that is a subset of brain mask where coupling should be computed
#' @param type "regression" or "pca"
#' @param ref reference modality when type="pca" or dependent modality when type="regression"
#' @param fwhm Full width in voxels for FWHM
#' @param thresh Threshold for trimming Gaussian kernel weights
#' @param radius radius of neighborhood in number of voxels. Can be specified instead of fwhm. Default is NULL and overrides fwhm when specified as non-NULL. 
#' @param reverse If TRUE, calculates both regressions if type="regression" and length(files)==2, otherwise ignored
#' @param verbose TRUE for updates on computation, else FALSE
#' @param retimg If TRUE, return list of estimated coupling maps as nifti objects
#' @param outDir Full path to directory where maps should be written
#' @param propMiss Maximum proportion of missing voxels in a neighborhood to tolerate, i.e., return NA if missing more than propMiss in the neighborhood of the center voxel
#' @export
#' @importFrom ANTsRCore antsGetSpacing getNeighborhoodInMask
#' @importFrom extrantsr check_ants 
#' @importFrom stats qnorm
#' @return Estimated IMCo coupling maps, either written to files and/or returned as nifti objects
imco <- function(files, brainMask, subMask=NULL, type="pca", ref=1, fwhm=3, 
                 thresh=0.005, radius=NULL, reverse=FALSE, 
                 verbose=TRUE, retimg=FALSE, outDir=NULL, propMiss=NULL){
    if(!(type=="pca" | type=="regression")){
        stop('type must be either pca or regression')
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
        if(!all(ANTsRCore::antsGetSpacing(fileList[[i-1]])==ANTsRCore::antsGetSpacing(fileList[[i]]))){
            stop('Voxel dimensions do not match')
        }
    }
    # Dimension of each voxel in mm
    vDims = ANTsRCore::antsGetSpacing(fileList[[1]])
    # Image dimension
    imgDims = dim(fileList[[1]])
    if(is.null(radius)){
    	if(fwhm <= 0){
    		stop('fwhm must be positive')
    	}
    	if(fwhm > min(imgDims)){
    		stop('fwhm too large')
    	}
    	if(!is.numeric(thresh) | thresh <= 0 | thresh >=1){
    		stop('thresh must be between 0 and 1')
    	} 
    } else{
    	if(radius < 1){
    		stop('radius must be 1 or greater')
    	}
    	if(radius > min(imgDims)){
    		stop('radius too large')
    	}
    	if(!is.numeric(thresh) | thresh <= 0 | thresh >=1){
    		stop('thresh must be between 0 and 1')
    	} 
    }
    # Read in brain mask
    bMask = extrantsr::check_ants(brainMask)
    if(!all(dim(bMask)==dim(fileList[[1]]))){
        stop('Image dimensions do not match the brain mask')
    }
    if(!all(ANTsRCore::antsGetSpacing(bMask)==ANTsRCore::antsGetSpacing(fileList[[1]]))){
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
    if(is.null(radius)){
    	minDim = min(vDims) #Should I use max or min to compute sigma?
    	width = fwhm*minDim
        # Need sigma for specifying Guassian kernel 
        # sigma is variance here
    	sigma = width/(2*sqrt(2*log(2)))
        # Backsolve for radius using kernel formula, assuming 3D vector with equal elements (e.g., c(2,2,2))
    	radium = round(sqrt(-2*sigma*log(thresh)/3), 0)
    	radius = rep(radium, 3)
    } else{
    	sigma = (radius^2)*(-3/2)/log(thresh)
    	radius = rep(radius, 3)
    }
    if(verbose){
    	cat("# Extracting neighborhood data \n")
    }
    # Neighborhood data from each modality
    if(!is.null(subMask)){
        sMask = extrantsr::check_ants(subMask)
        mask_indices = which(as.array(sMask) > 0)
        nhoods = lapply(fileList, function(x) ANTsRCore::getNeighborhoodInMask(image=x, mask=sMask, radius=radius, spatial.info=TRUE, boundary.condition='image'))
    } else{
        mask_indices = which(as.array(bMask) > 0)
        nhoods = lapply(fileList, function(x) ANTsRCore::getNeighborhoodInMask(image=x, mask=bMask, radius=radius, spatial.info=TRUE))
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
    	regObj = imco_reg(files=fileList, nhoods=nhoods, nWts=nWts, mask_indices=mask_indices, ref=ref, reverse=reverse, verbose=verbose, retimg=retimg, outDir=outDir, propMiss=propMiss)
    	return(regObj)
    } else{
    	pcaObj = imco_pca(files=fileList, nhoods=nhoods, nWts=nWts, mask_indices=mask_indices, ref=ref, verbose=verbose, retimg=retimg, outDir=outDir, propMiss=propMiss)
    	return(pcaObj)
    }
}


