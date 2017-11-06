#' @useDynLib bigtime, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @import Rcpp
#' @import methods
#' @import stats
#' @importFrom utils head
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom utils tail
#' @import corrplot
#' @import grDevices
#' @import graphics

.onUnload <- function (libpath) {
  library.dynam.unload("bigtime", libpath)
}


