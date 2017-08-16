#' @title make_ants_image
#'
#' @description Convert vector of IMCo into ANTS image
#' @param vec vector of IMCo values
#' @param mask_indices voxels where IMCo was computed
#' @param reference image for header info
#' @export
#' @return 3D ANTS image of IMCo measurements
#' @examples \dontrun{
#' 
#'}
make_ants_image = function(vec, mask_indices, reference){
	arr = array(0, dim=dim(reference))
	arr[mask_indices] = vec
	x = as.antsImage(arr, reference=reference)
	return(x)
}
