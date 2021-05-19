
#' Gives the residuals for VAR models estimated using sparseVAR
#'
#' @param object model estimated using sparseVAR
#' @param ... Not currently used
#' @export
#' @return Returns a matrix of residuals.
#' @examples
#' dat <- simVAR(200, 2, 5, decay = 0.001, seed = 6150533)
#' mod <- sparseVAR(dat$Y)
#' f <- fitted(mod)
#' res <- resid(mod)
residuals.bigtime.VAR <- function(object, ...){
  mod <- object
  fit <- fitted.bigtime.VAR(mod, ...)
  s <- dim(mod$Y)[1] - dim(fit)[1]
  res <- mod$Y[-(1:s), ] - fit
  colnames(res) <- colnames(mod$Y)
  res
}

#' Gives the residuals for VAR models estimated using sparseVAR
#'
#' @param object model estimated using sparseVAR
#' @param ... Not currently used
#' @export
#' @return Returns a matrix of residuals.
#' @examples
#' dat <- bigtime::simVAR(200, 2, 5, decay = 0.001, seed = 6150533)
#' mod <- bigtime::sparseVAR(dat$Y)
#' f <- fitted(mod)
#' res <- resid(mod)
resid.bigtime.VAR <- function(object, ...){
  residuals.bigtime.VAR(object, ...)
}


#' Gives the fitted values of a model estimated using sparseVAR
#'
#' @param object model estimated using sparseVAR
#' @param ... Not currently used
#' @export
#' @return Returns a matrix of fitted values
#' @examples
#' dat <- bigtime::simVAR(200, 2, 5, decay = 0.001, seed = 6150533)
#' mod <- bigtime::sparseVAR(dat$Y)
#' f <- fitted(mod)
#' res <- resid(mod)
#'
fitted.bigtime.VAR <- function(object, ...){
  mod <- object
  VARdata <- HVARmodel(Y=mod$Y, p=mod$p, h=mod$h)
  fit <- t(apply(VARdata$fullZ, 2, function(x) mod$Phihat%*%x + mod$phi0hat))
  colnames(fit) <- colnames(mod$Y)
  fit
}







