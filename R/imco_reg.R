#' @title Regression Inter-modal Coupling Model
#'
#' @description Implements full volumetric IMCo coupling estimation on a single subject using regression approach
#' @param files Vector of full paths to images or a list of local nifti objects
#' @param nhoods Neighborhood matrices computed in imco
#' @param nWts Vector of kernel weights
#' @param mask_indices Voxels where IMCo is computed
#' @param ref Which modality is dependent variable
#' @param reverse Should regression be done/returned both Y on X and X on Y? Ignored if # modalities > 2
#' @param verbose TRUE for updates on computation, else FALSE
#' @param retimg If TRUE, return list of estimated coupling maps as nifti objects
#' @param outDir Full path to directory where maps should be written
#' @export
#' @importFrom extrantsr ants2oro
#' @importFrom rlist list.rbind
#' @importFrom stats coef cov.wt lm.wfit
#' @return Estimated IMCo coupling maps, either written to files and/or returned as nifti objects
#' @examples \dontrun{
#'
#'}
imco_reg <- function(files, nhoods, nWts, mask_indices, ref=1, reverse=TRUE, verbose=TRUE, retimg=FALSE, outDir=NULL){
	if(verbose){
		cat("# Computing weighted regressions \n")
	}
	if(length(files)==2 & reverse){
		if(ref==1){
 			x = nhoods[[2]]$values
 			y = nhoods[[1]]$values
 		} else{
 			x = nhoods[[1]]$values
 			y = nhoods[[2]]$values
 		}
 		params = weighted_slr(x=x, y=y, wts=nWts, mInds=mask_indices, refImg=files[[1]])
 		if(!is.null(outDir)){
 			antsImageWrite(params$beta0, file.path(outDir, paste0('beta0.nii.gz')))
 			antsImageWrite(params$beta1, file.path(outDir, paste0('beta1.nii.gz')))
 			antsImageWrite(params$r2, file.path(outDir, paste0('rsquared.nii.gz')))
 			if(reverse){
 				antsImageWrite(params$beta0_reverse, file.path(outDir, paste0('beta0_reverse.nii.gz')))
 				antsImageWrite(params$beta1_reverse, file.path(outDir, paste0('beta1_reverse.nii.gz')))
 				antsImageWrite(params$r2_reverse, file.path(outDir, paste0('rsquared_reverse.nii.gz')))
 			}
 		}
 		if(retimg){
 			return(lapply(params, ants2oro))
 		} else{
 			return(NULL)
 		}
 	} else{
 		imgVals = lapply(nhoods, function(x) x$values)
 		bigDf = list.rbind(imgVals)
 		matList = lapply(split(bigDf, col(bigDf)), function(x) matrix(x, ncol=length(files)))
 		rmnaList = lapply(matList, function(x){
 			w=nWts
 			xRows = apply(as.matrix(x), 1, function(z){!any(is.na(z))})
 			if(sum(xRows) > 2){
 				return(cbind(w[xRows], as.matrix(x)[xRows,]))
 			}
 			return(NA)
 		})
 		wregList = lapply(rmnaList, function(x){
 			if(!is.na(x)[1]){
 				w=x[,1]
 				newx=x[,-1]
 				return(lm.wfit(x=as.matrix(cbind(1, newx[,-ref])),y=newx[,ref], w=w))
 			}
 			return(NA)
 		})
 		getR2 = lapply(wregList, function(x){
 			if(!is.na(x)[1]){
 				r = x$residuals
 				f = x$fitted
 				w = x$weights
 				wMean = sum(w * f/sum(w))
 				mss = sum(w * (f - wMean)^2)
 				rss = sum(w * (r^2))
 				r2 = mss/(mss + rss)
 				return(r2)
 			}
 			return(NA)
 		})
 		rsquared = make_ants_image(vec=as.vector(list.rbind(getR2)), mask_indices=mask_indices, reference=files[[1]])
 		if(!is.null(outDir)){
 			antsImageWrite(rsquared, file.path(outDir, 'rsquared.nii.gz'))
 		}
 		coefs = list()
 		for(j in 1:length(files)){
 			tmp = as.vector(list.rbind(lapply(wregList,
 				function(x){
 					if(!is.na(x)[1]) return(coef(x)[j]) else return(NA)
 				})
 			))
 			coefs[[j]] = make_ants_image(vec=tmp, mask_indices=mask_indices, reference=files[[1]])
 			if(!is.null(outDir)){
 				antsImageWrite(coefs[[j]], file.path(outDir, paste0('beta', j-1, '.nii.gz')))
 			}
 		}
 		if(retimg){
 			rsquared = ants2oro(rsquared)
 			intercept = ants2oro(coefs[[1]])
 			coefs[[1]] = NULL
 			slopes = lapply(coefs, ants2oro)
 			return(list('intercept'=intercept, 'slopes'=slopes, 'rsquared'=rsquared))
 		} else{
 			return(NULL)
 		}
 	}
}
