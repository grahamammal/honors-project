#' @useDynLib rsrHonors
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("rsrHonors", libpath)
}
