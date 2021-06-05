
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

  # We must catch the special case in which no selection precedure was used
  if (mod$selection == "none"){
    res <- array(dim = c(nrow(mod$Y)-mod$p, mod$k, length(mod$lambdas)))
    for (i in 1:length(mod$lambdas)){
      mod_tmp <- mod
      mod_tmp$Phihat <- mod$Phihat[, , i]
      mod_tmp$phi0hat <- mod$phi0hat[, , i]
      mod_tmp$selection <- "tmp"
      res[, , i] <- residuals.bigtime.VAR(mod_tmp, ...)
    }
    return(res)
  }

  fit <- fitted.bigtime.VAR(mod, ...)
  s <- dim(mod$Y)[1] - dim(fit)[1]
  res <- mod$Y[-(1:s), ] - fit
  colnames(res) <- colnames(mod$Y)
  res
}

#' Gives the residuals for VARX models estimated using sparseVARX
#'
#' @param object model estimated using sparseVARX
#' @param ... Not currently used
#' @export
#' @return Returns a matrix of residuals.
residuals.bigtime.VARX <- function(object, ...){
  mod <- object

  # Must cover the case when model was estimated using multiple lambda
  if (mod$selection == "none"){
    res <- array(NA, dim = c(nrow(mod$Y) - max(mod$p, mod$s), ncol(mod$Y), length(mod$lambdaPhi)))
    for (i in 1:length(mod$lambdaPhi)){
      mod_tmp <- mod
      mod_tmp$Phihat <- mod$Phihat[, , i]
      mod_tmp$phi0hat <- mod$phi0hat[, , i]
      mod_tmp$Bhat <- mod$Bhat[, , i]
      mod_tmp$selection <- "tmp"
      res[, , i] <- residuals.bigtime.VARX(mod_tmp, ...)
    }
    return(res)
  }
  fit <- fitted.bigtime.VARX(mod, ...)
  s <- nrow(mod$Y) - nrow(fit)
  res <- mod$Y[-(1:s), ] - fit
  colnames(res) <- colnames(mod$Y)
  res
}


#' Gives the residuals for VARMA models estimated using sparseVARMA
#'
#' @param object model estimated using sparseVARMA
#' @param ... Not currently used
#' @export
#' @return Returns a matrix of residuals.
residuals.bigtime.VARMA <- function(object, ...){
  mod <- object
  fit <- fitted.bigtime.VARMA(mod, ...)
  s <- nrow(mod$Y) - nrow(fit)
  res <- mod$Y[-(1:s), ] - fit
  colnames(res) <- colnames(mod$Y)
  res
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

  # We must catch the special case in which no selection procedure was used
  if (mod$selection == "none"){
    fit <- array(dim = c(nrow(mod$Y)-mod$p, mod$k, length(mod$lambdas)))
    for (i in 1:length(mod$lambdas)){
      mod_tmp <- mod
      mod_tmp$Phihat <- mod$Phihat[, , i]
      mod_tmp$phi0hat <- mod$phi0hat[, , i]
      mod_tmp$selection <- "tmp"
      fit[, , i] <- fitted.bigtime.VAR(mod_tmp, ...)
    }
    return(fit)
  }

  VARdata <- HVARmodel(Y=mod$Y, p=mod$p, h=mod$h)
  fit <- t(apply(VARdata$fullZ, 2, function(x) mod$Phihat%*%x + mod$phi0hat))
  if (ncol(mod$Y) == 1) fit <- t(fit) # somehow in this case we do not actually need to transpose
  colnames(fit) <- colnames(mod$Y)
  fit
}

#' Gives the fitted values of a model estimated using sparseVARX
#'
#' @param object model estimated using sparseVARX
#' @param ... Not currently used
#' @export
#' @return Returns a matrix of fitted values
fitted.bigtime.VARX <- function(object, ...){
  mod <- object

  # Must cover the case when model was estimated using multiple lambda
  if (mod$selection == "none"){
    fit <- array(NA, dim = c(nrow(mod$Y) - max(mod$p, mod$s), ncol(mod$Y), length(mod$lambdaPhi)))
    for (i in 1:length(mod$lambdaPhi)){
      mod_tmp <- mod
      mod_tmp$Phihat <- mod$Phihat[, , i]
      mod_tmp$phi0hat <- mod$phi0hat[, , i]
      mod_tmp$Bhat <- mod$Bhat[, , i]
      mod_tmp$selection <- "tmp"
      fit[, , i] <- fitted.bigtime.VARX(mod_tmp, ...)
    }
    return(fit)
  }

  dat <- HVARXmodel(mod$Y, mod$X, mod$p, mod$s, mod$h)
  X <- t(dat$fullX)
  Y <- t(dat$fullZ)

  fit <- lapply(1:nrow(X), function(i) mod$Phihat%*%Y[i, ] + mod$Bhat%*%X[i, ] + mod$phi0hat)
  fit <- t(do.call(cbind, fit))
  if (!is.null(colnames(mod$Y))) colnames(fit) <- colnames(mod$Y)
  fit
}


#' Gives the fitted values of a model estimated using sparseVARMA
#'
#' @param object model estimated using sparseVARMA
#' @param ... Not currently used
#' @export
#' @return Returns a matrix of fitted values
fitted.bigtime.VARMA <- function(object, ...){
  mod <- object

  # A VARMA is basically a VARX if we consider the errors as exogenous
  # We thus only need to restructure the model and can then use
  # fitted.bigtime.VARX

  mod_tmp <- mod
  mod_tmp$Bhat <- mod$Thetahat
  mod_tmp$X <- mod$U
  mod_tmp$p <- mod$VARMAp
  mod_tmp$s <- mod$VARMAq
  fitted.bigtime.VARX(mod_tmp)
}


#' Creates a diagnostic Plot
#'
#' @param mod VAR model estimated using sparseVAR
#' @param variable Variable to show. Either numeric (which column) or character
#' (variable name)
#' @param dates Optional Date vector.
#' @export
#' @return Returns a ggplot2 plot
#' @examples
#' dat <- simVAR(200, 2, 5, decay = 0.1, seed = 6150533,
#'                        sparsity_pattern = "hvar")
#' mod <- sparseVAR(dat$Y, selection = "bic", h = 1)
#' diagnostics_plot(mod, variable = 1)
diagnostics_plot <- function(mod, variable = 1, dates = NULL){
  UseMethod("diagnostics_plot", mod)
}


#' diagnostics_plot function for VAR models
#' @param mod VAR model estimated using sparseVAR
#' @param variable Variable to show. Either numeric (which column) or character
#' (variable name)
#' @param dates Optional Date vector.
#' @export
diagnostics_plot.bigtime.VAR <- function(mod, variable = 1, dates=NULL){
  if(!"bigtime.VAR" %in% class(mod)) stop("Only implemented for VAR models")
  if (mod$selection == "none") stop("Cannot be used with selection='none'. Choose other selection procedure in sparseVAR or call ic_selection on model first.")
  fit <- fitted.bigtime.VAR(mod)
  res <- residuals.bigtime.VAR(mod)
  s <- dim(mod$Y)[1] - dim(fit)[1]
  Y <- mod$Y[-(1:s), , drop = FALSE]

  .diagnostics_plot(Y, fit, res, s, variable, dates)
}


#' diagnostics_plot function for VARX models
#' @param mod VAR model estimated using sparseVARX
#' @param variable Variable to show. Either numeric (which column) or character
#' (variable name)
#' @param dates Optional Date vector.
#' @export
diagnostics_plot.bigtime.VARX <- function(mod, variable = 1, dates = NULL){
  if (!"bigtime.VARX" %in% class(mod)) stop("Only implemented for VARX models")
  if (mod$selection == "none") stop("Cannot be used with selection='none'. Choose other selection procedure in sparseVAR or call ic_selection on model first.")
  fit <- fitted.bigtime.VARX(mod)
  res <- residuals.bigtime.VARX(mod)
  s <- nrow(mod$Y) - nrow(fit)
  Y <- mod$Y[-c(1:s), , drop = FALSE]
  .diagnostics_plot(Y, fit, res, s, variable, dates)
}


#' diagnostics_plot function for VARMA models
#' @param mod VAR model estimated using sparseVARMA
#' @param variable Variable to show. Either numeric (which column) or character
#' (variable name)
#' @param dates Optional Date vector.
#' @export
diagnostics_plot.bigtime.VARMA <- function(mod, variable = 1, dates = NULL){
  if (!"bigtime.VARMA" %in% class(mod)) stop("Only implemented for VARMA models")
  fit <- fitted.bigtime.VARMA(mod)
  res <- residuals.bigtime.VARMA(mod)
  s <- nrow(mod$Y) - nrow(fit)
  Y <- mod$Y[-c(1:s), , drop = FALSE]
  .diagnostics_plot(Y, fit, res, s, variable, dates)
}


#' Internal function to plot the diagnostics plot
#' @param Y observed values
#' @param fit fitted values
#' @param res residual values
#' @param s how many observations were lost
#' @param variable variable to plot: either numeric or colname
#' @param dates vector of dates of length nrow(Y)
#' @keywords internal
.diagnostics_plot <- function(Y, fit, res, s, variable, dates){

  if (is.numeric(variable)) {
    if (variable > ncol(fit)) stop("Data does not have ", variable, " columns, but only ", ncol(fit))
    if (variable <= 0) stop("variable must be a name or an index >= 1")
  }else if(is.character(variable)) {
    if (is.null(colnames(fit))) stop("Data does not have column names.")
    else if (!(variable %in% colnames(fit))) stop("Data has no variable called ", variable)
    variable <- which(colnames(fit) == variable)
  }
  if (is.null(colnames(fit))){
    warning("Data has no column names. Using defaults Y1, ...")
    colnames(fit) <- paste0("Y", 1:ncol(fit))
    colnames(res) <- paste0("Y", 1:ncol(fit))
    colnames(Y) <-  paste0("Y", 1:ncol(fit))
  }
  variable <- colnames(fit)[variable]
  colnames(fit) <- paste0(colnames(fit), "_Fitted")
  colnames(Y) <- paste0(colnames(Y), "_Observed")
  colnames(res) <- paste0(colnames(res), "_Residual")
  df <- as.data.frame(cbind(fit, res, Y))

  d <- 1:nrow(fit) + s
  if(!is.null(dates)) {
    if (length(dates) != nrow(fit)) stop("dates must be of length ", nrow(fit), "not ", length(dates))
    d <- as.Date(dates)
  }

  cols <- c("Fitted"="red",
            "Observed"="black",
            "Residual"="black")
  Series <- Time <- Variable <- yintercept <- value <- p <- NULL # needed because of non-standard evaluation
  hline <- data.frame(list(p = "Residual", yintercept = 0))
  df %>%
    dplyr::mutate(Time = d) %>%
    tidyr::pivot_longer(-Time,
                        names_to = c("Variable", "Series"),
                        names_sep = "_") %>%
    dplyr::mutate(p = ifelse(Series %in% c("Fitted", "Observed"),
                             paste("Fitted", variable), "Residual"),
                  Series = factor(Series, levels = c("Observed",
                                                     "Fitted",
                                                     "Residual"))) %>%
    dplyr::filter(Variable == variable) %>%
    ggplot2::ggplot() +
    ggplot2::geom_hline(data = hline,
                        ggplot2::aes(yintercept = yintercept),
                        color = "red", linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(Time, value, color = Series)) +
    ggplot2::facet_grid(rows = dplyr::vars(p), scale = "free_y") +
    ggplot2::ylab("") +
    ggplot2::labs(color = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "top",
      legend.margin = ggplot2::margin(b=-10),
      # strip.background = ggplot2::element_blank(),
      # strip.text = ggplot2::element_blank()
    ) +
    ggplot2::scale_color_manual(values = cols, breaks = c("Observed", "Fitted"))
}


#' Checks whether a VAR is stable
#'
#' Using a model estimated by bigtime::sparseVAR, this function checks whether
#' the resulting VAR is stable
#'
#' @param mod model estimated using sparseVAR. Can only be a model
#' with one coefficient vector. Hence, the model must either be estimated using
#' cv=TRUE or by giving a single lambda value
#' @param verbose If TRUE, then the actual maximum eigenvalue of the companion
#' matrix will be printed to the console. Default is FALSE
#' @export
#' @return Returns TRUE if the VAR is stable and FALSE otherwise
is.stable <- function(mod, verbose = FALSE){
  if (!("bigtime.VAR") %in% class(mod)) stop("Model is not a VAR model estimated using sparseVAR")
  if (mod$selection == "none") stop("Model did not use any selection procedure. It is not clear which model is meant. Please use a selection procedure in sparseVAR or calls ic_selection on model.")
  Phi_hat <- mod$Phihat
  k <- mod$k
  p <- mod$p
  I <- diag(1, k*(p-1), k*(p-1))
  O <- matrix(0, k*(p-1), k)

  FF <- rbind(Phi_hat, cbind(I, O))
  max_eigval <- max(abs(eigen(FF)$value))
  if (verbose) cat("Maximum eigenvalue of Companion Matrix: ", max_eigval, "\n")
  if (max_eigval < 1) return(TRUE)
  FALSE
}


.check_if_standardised <- function(Y){
  Y_sd <- apply(Y, 2, sd)
  names(Y_sd) <- NULL
  Y_mu <- apply(Y, 2, mean)
  names(Y_mu) <- NULL
  if (!isTRUE(all.equal(rep(1, length(Y_sd)), Y_sd))) warning("It is recommended to standardise your data such that all variables have unit standard deviation")
  if (!isTRUE(all.equal(rep(0, length(Y_mu)), Y_mu))) warning("It is recommended to standardise your data such that all variables have zero mean.")
}








