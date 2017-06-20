res_dat <- function(.df) {
  tapply(.df[,4], .df[,1], length)
}
