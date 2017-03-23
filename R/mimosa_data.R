#' @title MIMoSA Training Data Frame
#'
#' @description This function creates the training vectors from a single MRI study that has FLAIR, T1, T2, and PD volumes as well as binary masks of lesions. The function can create a tissue mask for the data (or the user can supply a brain mask), the candidate voxels for lesion segmentation, smoothed volumes, and coupling maps. The user may supply already normalized data if they wish to use an alternative normalization method.
#' @param brain_mask brain or tissue mask of class nifti
#' @param FLAIR volume of class nifti
#' @param T1 volume of class nifti
#' @param T2 volume of class nifti. If not available use NULL.
#' @param PD volume of class nifti. If not available use NULL.
#' @param tissue is a logical value that determines whether the brain mask is a full brain mask or tissue mask (excludes CSF), should be FALSE unless you provide the tissue mask as the brain_mask object
#' @param normalize is a logical value that determines whether to perform z-score normalization of the image over the brain mask, should be TRUE unless you train model using an alternative normalization or provide normalized images
#' @param slices vector of desired slices to train on, if NULL then train over the entire brain mask
#' @param orientation string value telling which orientation the training slices are specified in, can take the values of "axial", "sagittal", or "coronal"
#' @param cores 1 numeric indicating the number of cores to be used (no more than 4 is useful for this software implementation)
#' @param verbose logical indicating printing diagnostic output
#' @export
#' @import fslr 
#' @import methods
#' @import nuerobase
#' @import oro.nifti
#' @import parallel
#' @import stats
#' @import oasis
#' @importFrom stats cov.wt qnorm
#' @return List of objects
#' @examples \dontrun{
#' 
#'}

mimosa_data <- function (brain_mask, FLAIR, T1, T2 = NULL, PD = NULL, tissue = FALSE, 
                                    gold_standard = NULL, normalize = TRUE, slices = NULL, 
                                    orientation = c("axial", "coronal", "sagittal"), cores = 1, verbose = TRUE) {
  if (verbose) {
    message("# Checking File inputs")
  }
  check_nifti2 <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    else {
      return(check_nifti(x))
    }
  }
  #verify nifti object
  FLAIR = check_nifti2(FLAIR)
  T1 = check_nifti2(T1)
  T2 = check_nifti2(T2)
  PD = check_nifti2(PD)
  gold_standard = check_nifti2(gold_standard)
  correct_image_dim2 <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    else {
      return(correct_image_dim(x))
    }
  }
  
  #correct image dim in case originals are weird
  FLAIR = correct_image_dim2(FLAIR)
  T1 = correct_image_dim2(T1)
  T2 = correct_image_dim2(T2)
  PD = correct_image_dim2(PD)
  gold_standard = correct_image_dim2(gold_standard)
  
  #if brain_mask is the tissue mask then just check nifti, make sure 0,1,and correct image dim
  #if brain_mask is not tissue mask then we create the tissue mask
  if(tissue == TRUE){
    brain_mask = check_nifti(brain_mask)
    brain_mask = brain_mask > 0
    brain_mask = datatyper(brain_mask, trybyte = TRUE)
    brain_mask = correct_image_dim(brain_mask)
    tissue_mask = brain_mask
  } else if (tissue == FALSE){
    if(verbose){
      message("# Getting Tissue Mask")
    }
    ero_brain_mask = fslerode(brain_mask, kopts = "-kernel box 5x5x5", 
                              retimg = TRUE)
    tissue_mask = voxel_selection(flair = FLAIR, brain_mask = ero_brain_mask, cutoff = 0.15)
  }
  #make a list of the images
  mimosa_data = list(FLAIR = FLAIR, T1 = T1, T2 = T2, 
                     PD = PD)
  rm(list = c("FLAIR","T1", "T2", "PD", "brain_mask"))
  
  #remove null values in case T2 or PD is null
  nulls = sapply(mimosa_data, is.null)
  mimosa_data = mimosa_data[!nulls]
  
  #normalize by z-score approach if specified
  if (normalize == TRUE) {
    if (verbose) {
      message("# Normalizing Images using Z-score")
    }
    mimosa_data = lapply(mimosa_data, zscore_img, mask = tissue_mask, 
                         margin = NULL)
    
    normalized = mimosa_data
    nulls = sapply(normalized, is.null)
    normalized = normalized[!nulls]
  }
  
  #obtain candidate voxels/candidate mask
  if (verbose) {
    message("# Voxel Selection Procedure")
  }
  top_voxels = voxel_selection(flair = mimosa_data$FLAIR, 
                               brain_mask = tissue_mask, cutoff = 0.85)
  if (verbose) {
    message("# Smoothing Images: Sigma = 10")
  }
  #obtain smoothed at 10 images
  smooth_10 = mclapply(mimosa_data, fslsmooth, sigma = 10, mask = tissue_mask, 
                       retimg = TRUE, smooth_mask = TRUE, mc.cores = cores)
  if (verbose) {
    message("# Smoothing Images: Sigma = 20")
  }
  #obtain smoothed at 20 images
  smooth_20 = mclapply(mimosa_data, fslsmooth, sigma = 20, mask = tissue_mask, 
                       retimg = TRUE, smooth_mask = TRUE, mc.cores = cores)
  #add smoothed images to the mimosa_data list with proper names
  img_names = names(mimosa_data)
  names(smooth_10) = paste0(img_names, "_10")
  mimosa_data = c(mimosa_data, smooth_10)
  names(smooth_20) = paste0(img_names, "_20")
  mimosa_data = c(mimosa_data, smooth_20)
  #remove the smoothed objects to conserve memory
  smoothed = list(smooth_10 = smooth_10, smooth_20 = smooth_20)
  rm(list = c("smooth_10", "smooth_20"))

  #initialize empty lists to store all the combinations of coupling intercepts and slopes
  coupling_intercepts = list()
  coupling_slopes = list()
  combos = combn(img_names, 2)
  combos = cbind(combos, rbind(combos[2,], combos[1,]))
  list_names = as.data.frame(matrix(nrow = 1, ncol = dim(combos)[2]))
  
  #run coupling on all possible combinations of images
  if (verbose) {
    message("# Running Coupling")
  }
  
  for(i in 1:dim(combos)[2]){
    
    temp_files = list(eval(parse(text = paste0('mimosa_data$', combos[1,i]))), 
                      eval(parse(text = paste0('mimosa_data$', combos[2,i]))))
    temp_return = imco(files=temp_files, brainMask=tissue_mask, subMask=top_voxels, type="regression", 
                       ref=1, neighborhoodSize=3, verbose=TRUE, retimg=TRUE, outDir=NULL)
    
    # Regress Y on X so var names will be YonX_int and YonX_slope
    list_names[i] = paste0(combos[1,i], 'on' , combos[2,i])
    coupling_intercepts[[i]] = temp_return$intercept
    coupling_slopes[[i]] = temp_return$slopes[[1]]
    
    if (verbose) {
      message(paste0('# Ran Coupling for ', combos[1,i], ' on ' , combos[2,i], ' Successfuly'))
    }
  }
  
  #rename coupling lists before appending with mimosa_data
  names(coupling_intercepts) = paste0(list_names, '_intercepts')
  names(coupling_slopes) = paste0(list_names, '_slopes')
  
  #append data
  mimosa_data = append(mimosa_data, coupling_intercepts)
  mimosa_data = append(mimosa_data, coupling_slopes)

  mimosa_data$gold_standard = gold_standard
  mimosa_data$top_voxels = top_voxels
  #turn full list into dataframe
  mimosa_dataframe = lapply(mimosa_data, function(X) masked = c(X[top_voxels == 1]))
  mimosa_dataframe = do.call(cbind, mimosa_dataframe)
  rownames(mimosa_dataframe) = NULL
  #now get indices of the top voxels to bind with the data for reference
  inds = niftiarr(top_voxels, 1)*top_voxels
  inds = which(inds > 0, arr.ind = TRUE)
  colnames(inds) = c("axial", "coronal", "sagittal")
  mimosa_dataframe = as.data.frame(cbind(inds, mimosa_dataframe))
  mimosa_dataframe$top_voxels = NULL
  
  #allow for slice images
  if (!is.null(slices)) {
    orientation = match.arg(orientation)
    mimosa_dataframe = mimosa_dataframe[mimosa_dataframe[, orientation] %in% slices, ]
  }

  #determine return object based on inputs
  if(normalize == FALSE & tissue == TRUE){
    #if normalize = FALSE we do not normalize, if tissue = true we treat the brain mask as the tissue_mask
    ##in this case do not return normalized images or the tissue mask as user already has them
    return(list(mimosa_dataframe = mimosa_dataframe, top_voxels = top_voxels, smoothed = smoothed, 
                coupling_intercepts = coupling_intercepts, coupling_slopes = coupling_slopes))
  }
  
  if(normalize == TRUE & tissue == TRUE){
    #if normalize is true then we normalize the images provided, if tissue is true we treat the brain mask as
    ##the tissue mask in this case return the normalized images but they have the tissue mask so do not return
    return(list(mimosa_dataframe = mimosa_dataframe, top_voxels = top_voxels, smoothed = smoothed, 
                coupling_intercepts = coupling_intercepts, coupling_slopes = coupling_slopes, 
                normalized = normalized))
  }
  if(normalize == FALSE & tissue == FALSE){
    #if normalize is FALSE then we normalize images if tissue is false we find the tissue mask
    ##return only tissue
    return(list(mimosa_dataframe = mimosa_dataframe, top_voxels = top_voxels, smoothed = smoothed, 
                coupling_intercepts = coupling_intercepts, coupling_slopes = coupling_slopes, 
                tissue_mask = tissue_mask))
    }
  if(normalize == TRUE & tissue == FALSE){
    #if normalize is true then images are normalized, if tissue is false then we find the tissue mask
    ##return both
    return(list(mimosa_dataframe = mimosa_dataframe, top_voxels = top_voxels, smoothed = smoothed, 
                coupling_intercepts = coupling_intercepts, coupling_slopes = coupling_slopes, 
                normalized = normalized, tissue_mask = tissue_mask))  
    }
}

  
  