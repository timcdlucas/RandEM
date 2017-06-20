bsD <- function(dat, tm, v){
  bsDat <- dat[sample(1:nrow(dat), size = nrow(dat), replace = TRUE), ]
  dOut <- rem(bsDat, tm, v)
  return(dOut)
}
