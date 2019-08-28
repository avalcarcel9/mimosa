#' @title Calculate performance statistics
#'
#' @description This function calculates true positive rate, false positive rate, false negative rate, false positive count, and sensitivity
#' @param gold_standard Gold standard segmentation mask of class nifti
#' @param predicted_segmentation Predicted segmentation mask volume of class nifti
#' @param k Minimum number of voxels for a segmentation cluster/component
#' @param percent_overlap Proportion of gold standard segmentation to be overlapped by predicted
#' @param verbose Logical indicating printing diagnostic output
#' @export
#' @importFrom extrantsr label_mask
#' @importFrom neurobase niftiarr
#' @return Matrix of results

count_stats <- function(gold_standard, predicted_segmentation, k = Inf, percent_overlap = NULL, verbose = TRUE){

  #Initialize results matrix
  stats_mat = matrix(NA, nrow = 1, ncol = 5)
  colnames(stats_mat) = c('True Positive Rate', 'False Positive Rate', 'False Negative Rate',
                          'False Positive Count', 'Sensitivity')

  #####################
  # Obtain Count Mask Images
  #####################
  gold_standard = check_nifti(gold_standard)
  predicted_segmentation = check_nifti(predicted_segmentation)

  message("# Counting Masks")

  gold_count = label_mask(gold_standard, k = k)
  predicted_count = label_mask(predicted_segmentation, k = k)

  #####################
  # True Positive
  #####################

  message("# Calculating True Positive")

  #Create Count Mask Dataframes
  inds = niftiarr(gold_standard, 1)*gold_count
  inds = which(inds > 0, arr.ind = TRUE)
  colnames(inds) = c("axial", "coronal", "sagittal")

  gold_count_df = lapply(list(gold_count = gold_count, predicted_segmentation = predicted_segmentation),
                         function(X) masked = c(X[gold_count > 0]))
  gold_count_df = cbind(inds, gold_count_df$gold_count, gold_count_df$predicted_segmentation)
  colnames(gold_count_df) = c("axial", "coronal", "sagittal", "gold_count", "predicted_indicator")

  gold_count_list = list()

  for(i in 1:max(gold_count_df[,4])){
  gold_count_list[[i]] = subset(gold_count_df, gold_count_df[,4]==i)
  }

  ############################
  # If percent_overlap is NULL then any overlapping voxel is considered true positive
  # If percent_overlap is proportion between 0 and 1 then lesion only is considered true positive if
  ## sufficient number of voxels are labelled lesion
  ############################

  if(is.null(percent_overlap)){
    TP <- function(count_mat){
      return(ifelse(sum(count_mat[,5])>0, 1, 0))
    }
  }
  else if(percent_overlap > 0 & percent_overlap <= 1){
    TP <- function(count_mat){
      return(ifelse(sum(count_mat[,5])>0 & sum(count_mat[,5])/dim(count_mat)[1] > percent_overlap, 1, 0))
    }
  }

  TP_vector = unlist(lapply(gold_count_list, TP))

  stats_mat[,1] = sum(TP_vector)/length(gold_count_list)

  #####################
  # False Negative
  #####################

  message("# Calculating False Negative")

  FN <- function(count_mat){
    return(ifelse(sum(count_mat[,5]) == 0, 1, 0))
  }

  FN_vector = unlist(lapply(gold_count_list, FN))

  stats_mat[,3] = sum(FN_vector)/length(gold_count_list)

  rm(inds, gold_count_df, gold_count_list, TP_vector)

  #####################
  # False Positive
  #####################

  message("# Calculating False Positive")

  #Create Count Mask Dataframes
  inds = niftiarr(predicted_segmentation, 1)*predicted_count
  inds = which(inds > 0, arr.ind = TRUE)
  colnames(inds) = c("axial", "coronal", "sagittal")

  pred_count_df = lapply(list(predicted_count = predicted_count, gold_standard = gold_standard),
                         function(X) masked = c(X[predicted_count > 0]))
  pred_count_df = cbind(inds, pred_count_df$predicted_count, pred_count_df$gold_standard)
  colnames(pred_count_df) = c("axial", "coronal", "sagittal", "pred_count", "GS_indicator")

  pred_count_list = list()

  for(i in 1:max(pred_count_df[,4])){
    pred_count_list[[i]] = subset(pred_count_df, pred_count_df[,4]==i)
  }

  #We miss the entire lesion (no voxels present in GS)
  FP <- function(count_mat){
    return(ifelse(sum(count_mat[,5]) == 0, 1, 0))
  }

  FP_vector = unlist(lapply(pred_count_list, FP))

  stats_mat[,2] = sum(FP_vector)/length(pred_count_list)

  stats_mat[,4] = sum(FP_vector)

  rm(inds, pred_count_df, pred_count_list, FP_vector)

  #####################
  # Sensitivity
  #####################

  message("# Calculating Sensitivity")

  stats_mat[,5] = stats_mat[,1] / (stats_mat[,1] + stats_mat[,3])

  return(as.data.frame(stats_mat))
}
