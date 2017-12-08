#' @title MIMoSA Training Data Frame
#'
#' @description This function creates the training vectors from a single MRI study that has FLAIR, T1, T2, and PD volumes as well as binary masks of lesions. The function can create a tissue mask for the data (or the user can supply a brain mask), the candidate voxels for lesion segmentation, smoothed volumes, and coupling maps. The user may supply already normalized data if they wish to use an alternative normalization method.
#' @param brain_mask brain or tissue mask of class nifti
#' @param FLAIR volume of class nifti
#' @param T1 volume of class nifti
#' @param T2 volume of class nifti. If not available use NULL.
#' @param PD volume of class nifti. If not available use NULL.
#' @param tissue is a logical value that determines whether the brain mask is a full brain mask or tissue mask (excludes CSF), should be FALSE unless you provide the tissue mask as the brain_mask object
#' @param gold_standard gold standard lesion segmentation mask of class nifti
#' @param normalize is NULL by default and will not perform any normalization on data. To normalize data specifcy Z for z-score normalization or WS for WhiteStripe normalization
#' @param cand_mask is NULL to use candidate mask procedure proposed with method or a nifti object to be used as the candidate mask
#' @param slices vector of desired slices to train on, if NULL then train over the entire brain mask
#' @param orientation string value telling which orientation the training slices are specified in, can take the values of "axial", "sagittal", or "coronal"
#' @param cores 1 numeric indicating the number of cores to be used (no more than 4 is useful for this software implementation)
#' @param verbose logical indicating printing diagnostic output
#' @export
#' @importFrom fslr fslerode fsl_smooth
#' @importFrom neurobase check_nifti datatyper niftiarr zscore_img
#' @importFrom parallel mclapply
#' @importFrom oasis voxel_selection correct_image_dim
#' @importFrom stats cov.wt qnorm
#' @importFrom utils combn
#' @importFrom WhiteStripe whitestripe whitestripe_norm
#' @return List of objects
#' @examples \dontrun{
#'
#'}

mimosa_data <- function (brain_mask, FLAIR, T1, T2 = NULL, PD = NULL, tissue = FALSE,
                                    gold_standard = NULL, normalize = FALSE, cand_mask = NULL, slices = NULL,
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
  # Verify inputs are NIFTI objects
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

  # Correct image dimension in case originals do not match
  FLAIR = correct_image_dim2(FLAIR)
  T1 = correct_image_dim2(T1)
  T2 = correct_image_dim2(T2)
  PD = correct_image_dim2(PD)
  gold_standard = correct_image_dim2(gold_standard)

  # Assignment or creation of tissue mask and verification mask is 0/1
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
  # Make a list of the images
  mimosa_data = list(FLAIR = FLAIR, T1 = T1, T2 = T2,
                     PD = PD)
  rm(list = c("FLAIR","T1", "T2", "PD", "brain_mask"))

  # Remove null values from list in case T2 or PD is specified as NULL
  nulls = sapply(mimosa_data, is.null)
  mimosa_data = mimosa_data[!nulls]

  # Normalize by z-score approach if normalize = "Z"
  if (normalize == 'Z') {
    if (verbose) {
      message("# Normalizing Images using Z-score")
    }
    mimosa_data = lapply(mimosa_data, zscore_img, mask = tissue_mask,
                         margin = NULL)

    normalized = mimosa_data
    nulls = sapply(normalized, is.null)
    normalized = normalized[!nulls]
  } 
  if (normalize == 'WS'){
    if (verbose) {
      message("# Normalizing Images using WhiteStripe")
    }
    # Run WhiteStripe on T1 and FLAIR
    WS_T1 = whitestripe(mimosa_data$T1, type="T1", stripped = TRUE, verbose = verbose)
    WS_T2 = whitestripe(mimosa_data$FLAIR, type="T2", stripped = TRUE, verose = verbose)
    
    mimosa_data$T1 = whitestripe_norm(mimosa_data$T1, WS_T1$whitestripe.ind)
    mimosa_data$FLAIR = whitestripe_norm(mimosa_data$FLAIR, WS_T2$whitestripe.ind)
    if(is.null(T2)==FALSE & is.null(PD)==TRUE){
      # Images include T1, FLAIR, T2
      mimosa_data$T2 = whitestripe_norm(mimosa_data$T2, WS_T2$whitestripe.ind)
    }
    if(is.null(T2)==TRUE & is.null(PD)==FALSE){
      # Images include T1, FLAIR, PD
      mimosa_data$PD = whitestripe_norm(mimosa_data$PD, WS_T2$whitestripe.ind)
    }
    if(is.null(T2)==FALSE & is.null(PD)==FALSE){
      # Images include T1, FLAIR, T2, PD
      mimosa_data$T2 = whitestripe_norm(mimosa_data$T2, WS_T2$whitestripe.ind)
      mimosa_data$PD = whitestripe_norm(mimosa_data$PD, WS_T2$whitestripe.ind)
    }
  }
  
  # Obtain candidate voxels and create candidate mask
  if (verbose) {
    message("# Voxel Selection Procedure")
  }
  if(is.null(cand_mask)){
    top_voxels = voxel_selection(flair = mimosa_data$FLAIR,
                                 brain_mask = tissue_mask, cutoff = 0.85)
  } else if(!is.null(cand_mask)){
    cand_mask = check_nifti(cand_mask)
    top_voxels = cand_mask
  }
  
  if (verbose) {
    message("# Smoothing Images: Sigma = 10")
  }
  # Obtain smoothed at 10 images
  
  smoothed_mask = fsl_smooth(file = tissue_mask, sigma = 10,
                                   smooth_mask = FALSE)
  
  smooth_10 = mclapply(mimosa_data, fslsmooth, sigma = 10, mask = tissue_mask,
                       retimg = TRUE, smooth_mask = TRUE, smoothed_mask = smoothed_mask, mc.cores = cores)
  if (verbose) {
    message("# Smoothing Images: Sigma = 20")
  }
  # Obtain smoothed at 20 images
  smoothed_mask = fsl_smooth(file = tissue_mask, sigma = 20,
                             smooth_mask = FALSE)
  
  smooth_20 = mclapply(mimosa_data, fslsmooth, sigma = 20, mask = tissue_mask,
                       retimg = TRUE, smooth_mask = TRUE, smoothed_mask = smoothed_mask, mc.cores = cores)
  
  # Add smoothed images to the mimosa_data list with proper names
  img_names = names(mimosa_data)
  names(smooth_10) = paste0(img_names, "_10")
  mimosa_data = c(mimosa_data, smooth_10)
  names(smooth_20) = paste0(img_names, "_20")
  mimosa_data = c(mimosa_data, smooth_20)
  # Remove the smoothed objects to conserve memory
  smoothed = list(smooth_10 = smooth_10, smooth_20 = smooth_20)
  rm(list = c("smooth_10", "smooth_20"))

  # Initialize empty lists to store all the combinations of coupling intercepts and slopes
  coupling_intercepts = list()
  coupling_slopes = list()
  combos = combn(img_names, 2)
  list_names = as.data.frame(matrix(nrow = 1, ncol = dim(combos)[2]*2))

  # Run coupling on all possible combinations of images
  if (verbose) {
    message("# Running Coupling")
  }

  for(i in 1:dim(combos)[2]){

    temp_files = list(eval(parse(text = paste0('mimosa_data$', combos[1,i]))),
                      eval(parse(text = paste0('mimosa_data$', combos[2,i]))))
    
    temp_return = imco(files=temp_files, brainMask=tissue_mask, subMask=top_voxels, type="regression",
                       ref=1, neighborhoodSize=3, reverse=TRUE, verbose=verbose, retimg=TRUE, outDir=NULL)

    # Regressed Y on X and X on Y so create variable names to match: YonX_int, YonX_slope, XonY_int, XonY_slope
    list_names[i] = paste0(combos[1,i], 'on' , combos[2,i])
    list_names[(i+dim(combos)[2])] = paste0(combos[2,i], 'on' , combos[1,i])
    coupling_intercepts[[i]] = temp_return$beta0
    coupling_slopes[[i]] = temp_return$beta1
    coupling_intercepts[[(i+dim(combos)[2])]] = temp_return$beta0_reverse
    coupling_slopes[[(i+dim(combos)[2])]] = temp_return$beta1_reverse

    if (verbose) {
      message(paste0('# Ran Coupling for ', combos[1,i], ' and ' , combos[2,i], ' Combinations Successfully'))
    }
  }

  # Rename coupling lists before appending with mimosa_data
  names(coupling_intercepts) = paste0(list_names, '_intercepts')
  names(coupling_slopes) = paste0(list_names, '_slopes')

  # Append data
  mimosa_data = append(mimosa_data, coupling_intercepts)
  mimosa_data = append(mimosa_data, coupling_slopes)

  mimosa_data$gold_standard = gold_standard
  mimosa_data$top_voxels = top_voxels
  # Transform full list of images into a dataframe
  mimosa_dataframe = lapply(mimosa_data, function(X) masked = c(X[top_voxels == 1]))
  mimosa_dataframe = do.call(cbind, mimosa_dataframe)
  rownames(mimosa_dataframe) = NULL
  # Bind indices of the top voxels to mimoda_data for reference
  inds = niftiarr(top_voxels, 1)*top_voxels
  inds = which(inds > 0, arr.ind = TRUE)
  colnames(inds) = c("axial", "coronal", "sagittal")
  mimosa_dataframe = as.data.frame(cbind(inds, mimosa_dataframe))
  mimosa_dataframe$top_voxels = NULL

  # Allow for slices of images
  if (!is.null(slices)) {
    orientation = match.arg(orientation)
    mimosa_dataframe = mimosa_dataframe[mimosa_dataframe[, orientation] %in% slices, ]
  }

  # Determine return object based on inputs
  if(normalize == FALSE & tissue == TRUE){
    # If normalize = FALSE we do not normalize, if tissue = true we treat the brain mask as the tissue_mask
    ## In this case do not return normalized images or the tissue mask as user already has them
    return(list(mimosa_dataframe = mimosa_dataframe, top_voxels = top_voxels, smoothed = smoothed,
                coupling_intercepts = coupling_intercepts, coupling_slopes = coupling_slopes))
  }

  if(normalize != FALSE & tissue == TRUE){
    # If normalize a value then we normalize the images provided, if tissue is true we treat the brain mask as
    ## The tissue mask in this case return the normalized images but they have the tissue mask so do not return
    return(list(mimosa_dataframe = mimosa_dataframe, top_voxels = top_voxels, smoothed = smoothed,
                coupling_intercepts = coupling_intercepts, coupling_slopes = coupling_slopes,
                normalized = normalized))
  }
  if(normalize == FALSE & tissue == FALSE){
    # If normalize is FALSE then we normalize images if tissue is false we find the tissue mask
    ## Return only tissue
    return(list(mimosa_dataframe = mimosa_dataframe, top_voxels = top_voxels, smoothed = smoothed,
                coupling_intercepts = coupling_intercepts, coupling_slopes = coupling_slopes,
                tissue_mask = tissue_mask))
    }
  if(normalize != FALSE & tissue == FALSE){
    # If normalize is true then images are normalized, if tissue is false then we find the tissue mask
    ## Return both
    return(list(mimosa_dataframe = mimosa_dataframe, top_voxels = top_voxels, smoothed = smoothed,
                coupling_intercepts = coupling_intercepts, coupling_slopes = coupling_slopes,
                normalized = normalized, tissue_mask = tissue_mask))
    }
}


