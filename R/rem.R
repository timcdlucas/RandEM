#' Random Encounter Model function
#' 
#' This function allows the automated calculation of animal densities.km2 using camera trap data and Random Encounter Modelling. 
#' 
#' The function assumes that the first 4 columns of the dataset contain: 
#' 1) An identifying number for each survey location (e.g. 1, 2, 3)
#' 2) The number of individuals of the focal species observed in each capture
#' 3) The radial distance to the detected animal in each capture, given in metres
#' 4) The angle of detection based on the location of the detected animal in each capture, given in radians
#' 
#' For a detailed example of how to use rem, click (\href{http://htmlpreview.github.io/?https://github.com/arcaravaggi/remBoot/blob/master/vignettes/remBoot.html}{here}).
#' 
#' @param dat = the data frame
#' @param tm (numeric) The total number of hours all cameras were left in-situ at a focal site
#' @param v (numeric) The distance travelled by the focal species in 24 hours, in kilometres
#' @keywords density, population, random encounter model
#' @export
#' @examples
#' 
#' data(hDat)
#' 
#' ## Split the data by survey site:
#' grpDat <- split_dat(hDat)
#' 
#' ## Define tm and v and pass the values to the function:
#' tm <- 3600
#' v <- 1.4
#' 
#' ## Use the rem function to calculate density estimates. For one survey site:
#' rem(hDat, tm, v)
#' ## Or
#' rem(hDat, tm = 3360, v = 1.4) 
#' 
#' ## For multiple survey sites, assuming tm and v are constant:
#' rem(dat = grpDat[[1]], tm, v)
#' 
#' ## If tm and v differ for each survey site, we can specify them alongside the REM function, as below. Note that if the focal species is a constant, v should not change.
#' 
#' rem(dat = grpDat[[1]], tm = 3600, v = 1.4)
#' rem(dat = grpDat[[2]], tm = 3360, v = 1.4) 
#' 
#' ## Before calculating variance, define the number of bootstrap iterations:
#' nboots <- 1000
#' 
#' ## Use the bootstrapping function boot_sd on each group dataframe n (i.e.nboots) times and calculate the standard deviation:
#' remsD <- lapply(grpDat, boot_sd) 
#' remsSD <- lapply(remsD, sd)
#' remsSD


rem <- function(dat, tm, v){
  # Calculate y: total captures
  y <- sum(dat[,2])
  # Calculate the first half of the equation above
  trm1 <- y/tm
  trm2 <- pi/((v*mean(dat[,3], 
                      na.rm = TRUE)) * (2 + (mean(dat[,4], 
                                                  na.rm = TRUE))))
  # Calculate D
  return(trm1 * trm2)
}
