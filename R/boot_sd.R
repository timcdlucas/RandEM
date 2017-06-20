#' Required for remsD bootstrapping
#' 
#' Resamples camera trap data to calculate variance for Random Encounter Model density estimates. Note that tm and v must be the same across all sites. If this is not the case, split your data before bootstrapping and run the function on each sub-group. Called by remsD. Calls bsD.
#' @param x = bsD resampled data
#' @export
#' @examples 
#' ## Define the number of bootstrapping iterations and apply boot_sd to the data:
#' tm <- 3600
#' v <- 0.85
#' nboots <- 1000
#' 
#' grpDat <- split_dat(hDat)
#' remsD <- lapply(grpDat, boot_sd)
#' remsSD <- lapply(remsD, sd)
#' remsSD
#' 
#' @importFrom stats quantile sd

boot_sd <-  function(x){
  d <- replicate(nboots, bsD(x, tm, v)) 
  return(d)
}
