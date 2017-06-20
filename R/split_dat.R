#' Data splitting for grpDat
#' 
#' Separate data into groups prior to bootstrapping. Calls grpDat_func.
#' 
#' @param .df = the data frame
#' @export
#' @examples
#' 
#' ## Read data:
#' data(hDat)
#' 
#' ## Split data by survey site:
#' grpDat <- split_dat(hDat)
#' 
#' ## Check that function has performed correctly:
#' lapply(grpDat, head)
#' lapply(grpDat, nrow)


split_dat <- function(.df) {
  lapply(unique(.df[[1]]), grpDat_func, df = .df)
}
