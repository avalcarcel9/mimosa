#' @title Train MIMoSA model on full training set
#'
#' @description This function trains the MIMoSA model from the data frames
#' produced by mimosa_data on all subjects and determines optimal threshold
#' based on training data
#' @param brain_mask vector of full path to brain mask
#' @param FLAIR vector of full path to FLAIR
#' @param T1 vector of full path to T1
#' @param T2 vector of full path to T2 if available. If not use NULL.
#' @param PD vector of full path to PD if available. If not use NULL.
#' @param tissue is a logical value that determines whether the brain mask is
#' a full brain mask or tissue mask (excludes CSF), should be FALSE unless
#' you provide the tissue mask as the brain_mask object
#' @param gold_standard vector of full path to Gold standard segmentations.
#' Typically manually segemented images.
#' @param normalize TRUE normalizes image inputs using z-score normalization
#' @param slices vector of desired slices to train on, if NULL then train
#' over the entire brain mask
#' @param orientation string value telling which orientation the training
#' slices are specified in, can take the values of "axial", "sagittal",
#' or "coronal"
#' @param cores numeric indicating the number of cores to be used
#' (no more than 4 is useful for this software implementation)
#' @param verbose logical indicating printing diagnostic output
#' @param outdir vector of paths/IDs to be pasted to objects that will be
#' saved. NULL if objects are not to be saved
#' @param optimal_threshold NULL. To run algorithm provide vector of thresholds
#' @export
#' @importFrom neurobase writenii readnii niftiarr
#' @importFrom data.table rbindlist
#' @importFrom stats predict
#' @importFrom fslr fslsmooth
#' @importFrom utils write.csv
#' @return GLM objects fit in the MIMoSA procedure and optimal threshold
#' evaluated for full training set
#' @examples \dontrun{
#'
#'}

mimosa_training <- function(brain_mask, FLAIR, T1, T2 = NULL, PD = NULL, tissue = FALSE,
                            gold_standard, normalize = TRUE, slices = NULL,
                            orientation = c("axial", "coronal", "sagittal"), cores = 1, verbose = TRUE, outdir = NULL,
                            optimal_threshold = NULL){

  if(!(all.equal(length(brain_mask), length(FLAIR), length(T1), length(T2), length(PD)))){
    stop('path vectors do not match')
  }

  if (verbose) {
    message('# Obtaining training subject data')
  }

  top_voxel_list = list()
  train_data_all_list = list()
  gs_DSC_list = list()

  for(i in 1:length(brain_mask)){

    subject_files = cbind(brain_mask[i], FLAIR[i], T1[i], T2[i], PD[i], gold_standard[i])

    #names/formula/mimosa_data for items in subject_files and formula to train based on T2/PD missingness
    if(is.null(T2)==TRUE & is.null(PD)==FALSE){
      # Images include T1, FLAIR, PD
      formula = gold_standard ~ FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + PD_10 * PD + PD_20 * PD + T1_10 * T1 + T1_20 * T1 +
        FLAIRonT1_intercepts + FLAIRonPD_intercepts + T1onPD_intercepts +
        T1onFLAIR_intercepts + PDonFLAIR_intercepts + PDonT1_intercepts +
        FLAIRonT1_slopes + FLAIRonPD_slopes + T1onPD_slopes +
        T1onFLAIR_slopes + PDonFLAIR_slopes + PDonT1_slopes
      if (all(file.exists(subject_files))) {
        imgs_list  = lapply(subject_files, readnii, reorient = FALSE)
      }
      names(imgs_list) = c('brain_mask', 'FLAIR', 'T1', 'PD', 'gold_standard')
      train_data_i = mimosa_data(brain_mask = imgs_list$brain_mask, FLAIR = imgs_list$FLAIR, T1 = imgs_list$T1, T2 = NULL, PD = imgs_list$PD,
                                 tissue = tissue, gold_standard = imgs_list$gold_standard, normalize = normalize, slices = slices,
                                 orientation = orientation, cores = cores, verbose = verbose)
    }
    if(is.null(T2)==FALSE & is.null(PD)==TRUE){
      # Images include T1, FLAIR, T2
      formula = gold_standard ~ FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + T2_10 * T2 + T2_20 * T2 + T1_10 * T1 + T1_20 * T1 +
        FLAIRonT1_intercepts + FLAIRonT2_intercepts + T1onT2_intercepts +
        T1onFLAIR_intercepts + T2onFLAIR_intercepts + T2onT1_intercepts +
        FLAIRonT1_slopes + FLAIRonT2_slopes + T1onT2_slopes +
        T1onFLAIR_slopes + T2onFLAIR_slopes + T2onT1_slopes
      if (all(file.exists(subject_files))) {
        imgs_list  = lapply(subject_files, readnii, reorient = FALSE)
      }
      names(imgs_list) = c('brain_mask', 'FLAIR', 'T1', 'T2', 'gold_standard')
      train_data_i = mimosa_data(brain_mask = imgs_list$brain_mask, FLAIR = imgs_list$FLAIR, T1 = imgs_list$T1, T2 = imgs_list$T2, PD = NULL,
                                 tissue = tissue, gold_standard = imgs_list$gold_standard, normalize = normalize, slices = slices,
                                 orientation = orientation, cores = cores, verbose = verbose)
    }
    if(is.null(T2)==TRUE & is.null(PD)==TRUE){
      # Images include T1 and FLAIR
      formula = gold_standard ~ FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + T1_10 * T1 + T1_20 * T1 +
        FLAIRonT1_intercepts + T1onFLAIR_intercepts +
        FLAIRonT1_slopes + T1onFLAIR_slopes
      if (all(file.exists(subject_files))) {
        imgs_list  = lapply(subject_files, readnii, reorient = FALSE)
      }
      names(imgs_list) = c('brain_mask', 'FLAIR', 'T1', 'gold_standard')
      train_data_i = mimosa_data(brain_mask = imgs_list$brain_mask, FLAIR = imgs_list$FLAIR, T1 = imgs_list$T1, T2 = NULL, PD = NULL,
                                 tissue = tissue, gold_standard = imgs_list$gold_standard, normalize = normalize, slices = slices,
                                 orientation = orientation, cores = cores, verbose = verbose)
    }
    if(is.null(T2)==FALSE & is.null(PD)==FALSE){
      # Images include T1, FLAIR, T2, PD
      formula = gold_standard ~ FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + PD_10 * PD + PD_20 * PD + T2_10 * T2 + T2_20 * T2 + T1_10 * T1 + T1_20 * T1 +
        FLAIRonT1_intercepts + FLAIRonT2_intercepts + FLAIRonPD_intercepts +
        T1onT2_intercepts + T1onPD_intercepts + T2onPD_intercepts +
        T1onFLAIR_intercepts + T2onFLAIR_intercepts + PDonFLAIR_intercepts +
        T2onT1_intercepts + PDonT1_intercepts + PDonT2_intercepts +
        FLAIRonT1_slopes + FLAIRonT2_slopes + FLAIRonPD_slopes +
        T1onT2_slopes + T1onPD_slopes + T2onPD_slopes +
        T1onFLAIR_slopes + T2onFLAIR_slopes + PDonFLAIR_slopes +
        T2onT1_slopes + PDonT1_slopes + PDonT2_slopes
      if (all(file.exists(subject_files))) {
        imgs_list  = lapply(subject_files, readnii, reorient = FALSE)
      }
      names(imgs_list) = c('brain_mask', 'FLAIR', 'T1', 'T2', 'PD', 'gold_standard')
      train_data_i = mimosa_data(brain_mask = imgs_list$brain_mask, FLAIR = imgs_list$FLAIR, T1 = imgs_list$T1, T2 = imgs_list$T2, PD = imgs_list$PD,
                                 tissue = tissue, gold_standard = imgs_list$gold_standard, normalize = normalize, slices = slices,
                                 orientation = orientation, cores = cores, verbose = verbose)
    }

    # Save candidate mask and mimosa_dataframe since we need it later in lists
    top_voxel_list[[i]] = train_data_i$top_voxels
    train_data_all_list[[i]] = train_data_i$mimosa_dataframe
    gs_DSC_list[[i]] = sum(imgs_list$gold_standard)

    if (verbose) {
      message(paste0('# Training Information for', brain_mask[i], 'Complete'))
    }

    if(!is.null(outdir)){

      if (verbose) {
        message(paste0('# Saving Subject Information for', outdir[i]))
      }

      # Write training dataframe for subject i
      ### Put outdir first so that we have the path as specified and ID

      # Return the train dataframe
      write.csv(train_data_i$mimosa_dataframe, file = paste0(outdir[i], '_mimosa_dataframe.csv'))
      # Write top voxels for subject i
      writenii(train_data_i$top_voxels, filename = paste0(outdir[i], '_top_voxels'))
      # Return the smoothed at 10 images
      for(j in 1:length(train_data_i$smoothed$smooth_10)){
        writenii(train_data_i$smoothed$smooth_10[[j]], filename = paste0(outdir[i], '_', names(train_data_i$smoothed$smooth_10)[j], '_smoothed'))
      }
      # Return the smoothed at 20 images
      for(j in 1:length(train_data_i$smoothed$smooth_20)){
        writenii(train_data_i$smoothed$smooth_20[[j]], filename = paste0(outdir[i], '_', names(train_data_i$smoothed$smooth_20)[j], '_smoothed'))
      }
      # Return the coupling intercept images
      for(j in 1:length(train_data_i$coupling_intercepts)){
        writenii(train_data_i$coupling_intercepts[[j]], filename = paste0(outdir[i], '_coupling_', names(train_data_i$coupling_intercepts)[j]))
      }
      # Return the slope images
      for(j in 1:length(train_data_i$coupling_slopes)){
        writenii(train_data_i$coupling_slopes[[j]], filename = paste0(outdir[i], '_coupling_', names(train_data_i$coupling_slopes)[j]))
      }
      ## Return normalized and/or tissue depending on inputs
      if(normalize == TRUE & tissue == TRUE){
        # If normalize is true then we normalize the images provided, if tissue is true we treat the brain mask as
        ##the tissue mask in this case return the normalized images but they have the tissue mask so do not return
        # Normalized images
        for(j in 1:length(train_data_i$normalized)){
          writenii(train_data_i$normalized[[j]], filename = paste0(outdir[i], '_', names(train_data_i$normalized)[j], '_norm'))
        }
      }
      if(normalize == FALSE & tissue == FALSE){
        # If normalize is FALSE then we normalize images if tissue is false we find the tissue mask
        ## Return tissue
        writenii(train_data_i$tissue_mask, filename = paste0(outdir[i], '_tissue_mask'))
      }
      if(normalize == TRUE & tissue == FALSE){
        # If normalize is true then images are normalized, if tissue is false then we find the tissue mask
        ## return both
        # Normalized images
        for(j in 1:length(train_data_i$normalized)){
          writenii(train_data_i$normalized[[j]], filename = paste0(outdir[i], '_', names(train_data_i$normalized)[j], '_norm'))
        }
        #tissue mask
        writenii(train_data_i$tissue_mask, filename = paste0(outdir[i], '_tissue_mask'))

      }
    }
  }

  # Transform list to dataframe so that we can fit the MIMoSA model
  train_dataframe_all = rbindlist(train_data_all_list)

  if (verbose) {
    message(paste0('# Fitting MIMoSA Model'))
  }

  # Fit Full MIMoSA Model
  mimosa_fit_model = mimosa_fit(training_dataframe=train_dataframe_all, formula=formula)
  rm(train_dataframe_all)

  # If user does not want to calculate optimal threshold (optimal_threshold == FALSE) then return the model
  if(is.null(optimal_threshold)){
    return(mimosa_fit_model)
  }
  if (!is.null(optimal_threshold)){
    
    # Initialize a storage matrix for DSC values
    dsc_mat = matrix(NA, nrow = length(train_data_all_list), ncol = length(optimal_threshold))

    for(i in 1:length(train_data_all_list)){

      # First generate probability maps for each subject
      predictions = predict(mimosa_fit_model, train_data_all_list[[i]], type = "response")
      predictions_nifti = niftiarr(top_voxel_list[[i]], 0)
      predictions_nifti[top_voxel_list[[i]] == 1] = predictions
      prob_map = fslsmooth(predictions_nifti, sigma = 1.25, mask = top_voxel_list[[i]],
                           retimg = TRUE, smooth_mask = TRUE)

      # Loop Through Thresholds to determine threshold
      DSC_scores = numeric()
      for (j in 1:length(optimal_threshold)) {
        # Threshold To Create Segmented Maps For Train Subjects
        lesion_mask = (prob_map >= optimal_threshold[j])
        lesion_df=c(lesion_mask[top_voxel_list[[i]] == 1])

        # New threshold1#
        DSC_scores[j] = (2*sum(lesion_df*train_data_all_list[[i]]$gold_standard))/(sum(lesion_df)+gs_DSC_list[[i]])
      }
      dsc_mat[i,] = DSC_scores
    }
    est_optimal_threshold=optimal_threshold[which.max(apply(dsc_mat, 2,mean))]
    return(list(mimosa_fit_model = mimosa_fit_model, estimated_optimal_threshold = est_optimal_threshold))

  }

}

