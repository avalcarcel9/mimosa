#' @title PCA Inter-modal Coupling Model
#'
#' @description Implements full volumetric IMCo coupling estimation on a single subject using PCA approach
#' @param files Vector of full paths to images or a list of local nifti objects
#' @param nhoods Neighborhood matrices computed in imco
#' @param nWts Vector of kernel weights
#' @param mask_indices Voxels where IMCo is computed
#' @param ref Which modality is dependent variable
#' @param verbose TRUE for updates on computation, else FALSE
#' @param retimg If TRUE, return list of estimated coupling maps as nifti objects
#' @param outDir Full path to directory where maps should be written
#' @param propMiss Maximum proportion of missing voxels in a neighborhood to tolerate, i.e., return NA if missing more than propMiss in the neighborhood of the center voxel
#' @export
#' @importFrom ANTsRCore antsImageWrite
#' @importFrom extrantsr ants2oro
#' @importFrom rlist list.rbind
#' @importFrom stats cov.wt
#' @return Estimated IMCo coupling maps, either written to files and/or returned as nifti objects
#' @examples \dontrun{
#' 
#'}
imco_pca <- function(files, nhoods, nWts, mask_indices, ref=1, verbose=TRUE, retimg=FALSE, outDir=NULL, propMiss=NULL){
    # Restructure to get eigen decomp at each voxel
	imgVals = lapply(nhoods, function(x) x$values)
	bigDf = list.rbind(imgVals)
	matList = lapply(split(bigDf, col(bigDf)), function(x) matrix(x, ncol=length(files)))
	rmnaList = lapply(matList, function(x){
		w=nWts
		xRows = apply(as.matrix(x), 1, function(z){!any(is.na(z))})
 			if(sum(xRows) > 2){
                # If the proportion of missing voxels is greater than propMiss, return NA
                if(mean(!xRows) > propMiss){
                   return(NA)
                } else{
                    return(cbind(w[xRows], as.matrix(x)[xRows,]))
                }
 			} 
		return(NA)
	})
	if(verbose){
		cat("# Computing weighted covariances \n")
	}
	rmnaListCenter = lapply(rmnaList, function(x){
		if(!is.na(x)[1]){
			w = x[,1]
			newx = scale(x[,-1], center=TRUE, scale=FALSE)
			return(cbind(w, newx))
		}
		return(NA)
	})
	rm(rmnaList)
    # Weighted cov of each matrix in matList
	wcovList = lapply(rmnaListCenter, function(x){
		if(!is.na(x)[1]){
			w=x[,1]
			newx = x[,-1]
			newx = cbind(newx[,ref], newx[,-ref])
			return(cov.wt(newx, wt=w, center=FALSE)$cov)
		}
		return(NA)
	})
	if(verbose){
		cat("# Computing weighted PCAs \n")
	}
	eigenList = lapply(wcovList, function(x){
		if(!is.na(x)[1]) return(eigen(x)) else return(NA)
	})
	rm(wcovList)
	if(verbose){
		cat("# Extracting IMCo images \n")
	}
	evals = list()
	components = list()
	for(j in 1:length(files)){
		tmp = as.vector(list.rbind(lapply(eigenList,
                function(x){
                    if(!is.na(x)[1]) return(x$values[j]) else return(NA)
                })
            ))
		evals[[j]] = make_ants_image(vec=tmp, mask_indices=mask_indices, reference=files[[1]])
		if(!is.null(outDir)){
			ANTsRCore::antsImageWrite(evals[[j]], file.path(outDir, paste0('eigenValue-', j, '.nii.gz')))    	
		}
		components[[j]] = list()
		for(k in 1:length(files)){
			tmp = as.vector(list.rbind(mapply(function(x, y){ 
				if(!is.na(x)[1]) return(x$vectors[k,j]) else return(NA)
			}, x=eigenList, y=rmnaListCenter)
			))
			components[[j]][[k]] = make_ants_image(vec=tmp, mask_indices=mask_indices, reference=files[[1]])
			if(!is.null(outDir)){
				ANTsRCore::antsImageWrite(components[[j]][[k]], file.path(outDir, paste0('component', j, '-', k, '.nii.gz')))    	
			}
		}
	}
	rm(eigenList)
	if(retimg){
		evals = lapply(evals, ants2oro)
		for(j in 1:length(files)){
			temp = lapply(components[[j]], ants2oro)
			components[[j]] = lapply(components[[j]], ants2oro)
		}
		return(list('eigenValueImages'=evals, 'eigenVectorImages'=components))
	} else{
		return(NULL)
	}
}
