#' @title Dataframe to Train MIMoSA Model
#'
#' @description Makes a dataframe to train, the tissue mask if desired, and coupling maps (intercepts and slopes) for all combos
#' @param brain_mask local brain mask object
#' @param FLAIR local FLAIR object
#' @param T1 local FLAIR object
#' @param T2 local T2 object if available. If not use NULL.
#' @param PD local PD object if available. If not use NULL.
#' @param tissue FALSE if brain mask is not the tissue mask (excludes CSF). TRUE if the brain_mask is the tissue mask.
#' @param gold_standard Gold standard segmentations. Typically manually segemented images.
#' @param normalize TRUE normalizes image inputs using z-score normalization
#' @param slices = NULL if full brain images are used
#' @param orientation c("axial", "coronal", "sagittal") orientation of slices
#' @param cores 1 is number of cores to be used
#' @param verbose TRUE allows progress update as function runs
#' @export
#' @import fslr 
#' @import methods
#' @import nuerobase
#' @import oro.nifti
#' @import parallel
#' @import stats
#' @import oasis
#' @importFrom stats cov.wt qnorm
#' @return List with dataframe prepared to be trained
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

  mimosa_data$GoldStandard = gold_standard
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
  cn = colnames(mimosa_dataframe)
  cn = setdiff(cn, orientation)
  mimosa_dataframe = mimosa_dataframe[, cn]
  #set return object
  if(tissue){
    return(list(mimosa_dataframe = mimosa_dataframe, voxel_selection = top_voxels, 
                gold_standard = mimosa_data$GoldStandard, coupling_intercepts = coupling_intercepts, 
                coupling_slopes = coupling_slopes))
  } else {
    return(list(mimosa_dataframe = mimosa_dataframe, voxel_selection = top_voxels, 
                gold_standard = mimosa_data$GoldStandard, coupling_intercepts = coupling_intercepts, 
                coupling_slopes = coupling_slopes, tissue_mask = tissue_mask))
  }
}