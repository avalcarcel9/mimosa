#' @title Inter-modal Coupling Model
#'
#' @description Implements full volumetric IMCo coupling estimation on a single subject
#' @param files Vector of full paths to images
#' @param brainMask Full path to brain mask
#' @param subMask NULL or full path subset of brain mask where coupling should be computed
#' @param type "regression" or "pca"
#' @param ref reference modality or dependent modality when type="regression"
#' @param neighborhoodSize Full width in voxels for FWHM
#' @param verbose TRUE for progress bar
#' @param retimg If TRUE, return list of estimated coupling maps
#' @param outDir Full path to directory where maps should be written
#' @export
#' @import ANTsR 
#' @import extrantsr 
#' @import rlist
#' @importFrom stats cov.wt qnorm
#' @return Estimated IMCo coupling maps, either written to files and/or returned as antsImage objects
#' @examples \dontrun{
#' 
#'}
imco <- function(files, brainMask, subMask=NULL, type="pca", ref=1, neighborhoodSize=3, verbose=TRUE, retimg=FALSE, outDir=NULL){
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
  ################################################
  # Restructure to get eigen decomp at each voxel
  ################################################
  imgVals = lapply(nhoods, function(x) x$values)
  bigDf = list.rbind(imgVals)
  matList = lapply(split(bigDf, col(bigDf)), function(x) matrix(x, ncol=length(files)))
  rmnaList = lapply(matList, function(x){w = nWts; xRows = apply(as.matrix(x), 1, function(z){!any(is.na(z))}); return(cbind(w[xRows], as.matrix(x)[xRows,]))})
  if(type=="regression"){
    if(verbose){
      cat("# Computing weighted regressions \n")
    }
    wregList = lapply(rmnaList, function(x){w=x[,1]; newx=x[,-1]; return(lm.wfit(x=as.matrix(cbind(1, newx[,-ref])),y=newx[,ref], w=w))})
    coefs = list()
    for(j in 1:length(files)){
      coefs[[j]] = array(0, dim=dim(fileList[[j]]))
      coefs[[j]][inds1] = as.vector(list.rbind(lapply(wregList, function(x) coef(x)[j])))
      coefs[[j]] = as.antsImage(coefs[[j]])
      coefs[[j]] = antsCopyImageInfo(target=coefs[[j]], reference=fileList[[j]])
      if(!is.null(outDir)){
        antsImageWrite(coefs[[j]], file.path(outDir, paste0('beta', j-1, '.nii.gz')))      
      }
    }
    if(retimg){
      intercept = ants2oro(coefs[[1]])
      coefs[[1]] = NULL
      slopes = lapply(coefs, ants2oro)
      return(list('intercept'=intercept, 'slopes'=slopes))
    }
  } else{
    if(verbose){
      cat("# Computing weighted covariances \n")
    }
    # Weighted cov of each matrix in matList
    wcovList = lapply(rmnaList, function(x){w=x[,1]; newx = x[,-1]; newx = cbind(newx[,ref], newx[,-ref]); return(cov.wt(newx, wt=w, center=TRUE)$cov)})
    if(verbose){
      cat("# Computing weighted PCAs \n")
    }
    eigenList = lapply(wcovList, function(x){return(eigen(x))})
    if(verbose){
      cat("# Extracting IMCo images \n")
    }
    evals = list()
    components = list()
    for(j in 1:length(files)){
      evals[[j]] = array(0, dim=dim(fileList[[j]]))
      evals[[j]][inds1] = as.vector(list.rbind(lapply(eigenList, function(x) x$values[j])))
      evals[[j]] = as.antsImage(evals[[j]])
      evals[[j]] = antsCopyImageInfo(target=evals[[j]], reference=fileList[[j]])
      if(!is.null(outDir)){
        antsImageWrite(evals[[j]], file.path(outDir, paste0('eigenValue-', j, '.nii.gz')))    	
      }
      components[[j]] = list()
      for(k in 1:length(files)){
        components[[j]][[k]] = array(0, dim=dim(fileList[[j]]))
        components[[j]][[k]][inds1] = as.vector(list.rbind(lapply(eigenList, function(x) x$vectors[k,j])))
        components[[j]][[k]] = as.antsImage(components[[j]][[k]])
        components[[j]][[k]] = antsCopyImageInfo(target=components[[j]][[k]], reference=fileList[[j]])
        if(!is.null(outDir)){
          antsImageWrite(components[[j]][[k]], file.path(outDir, paste0('component', j, '-', k, '.nii.gz')))    	
        }
      }
    }
    if(retimg){
      evals = lapply(evals, ants2oro)
      for(j in 1:length(files)){
        temp = lapply(components[[j]], ants2oro)
        components[[j]] = lapply(components[[j]], ants2oro)
      }
      return(list('eigenValueImages'=evals, 'eigenVectorImages'=components))
    }
  }
}