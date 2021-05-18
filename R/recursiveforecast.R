


#' Recursively Forecasts a VAR
#'
#' Recursively forecasts a VAR estimated using bigtime::sparseVAR.
#' lambda can either be NULL, in which case all lambdas that were used for
#' model estimation are used for forecasting, or a single value, in which case
#' only the model using this lambda will be used for forecasting.
#'
#' @param mod VAR model estimated using bigtime::sparseVAR
#' @param h forecast horizon: default=1
#' @param lambda Either NULL in which case a forecast will be made for all
#' lambdas for which the model was estimated, or a single value in which
#' case a forecast will only be made for the model using this lambda
#' @export
#' @return Returns an object of S3 class bigtime.recursiveforecast containing
#' \item{fcst}{Matrix or 3D array of forecasts}
#' \item{h}{Selected forecast horizon}
#' \item{lambda}{list of lambdas for which the forecasts were made}
#' \item{Y}{Data used for recursive forecasting}
recursiveforecast <- function(mod, h=1, lambda = NULL){
  if (!("bigtime.VAR" %in% class(mod))) stop("Recursive Forecasting is only supported for VAR models estimated using bigtime::sparseVAR")
  if (!is.null(lambda) & length(lambda) > 1) stop("Must either forecast using all lambdas or only one.")

  # If CV was not run, then PHI will be a cube.
  # in this case we can either use one slice or forecast for multiple lambda
  Phi_hat <- mod$Phihat
  phi_0 <- mod$phi0hat
  if (!is.matrix(Phi_hat) & is.null(lambda)) message("Forecasting using all lambdas")
  if (!is.null(lambda)){
    if (!is.matrix(Phi_hat) & !(lambda %in% mod$lambdas)) stop("Lambda must have been used for estimation.")
    if (is.matrix(Phi_hat)) warning("Choice of lambda will be ignored since model uses only the optimal lambda.")
  }
  if (!is.matrix(Phi_hat) & !is.null(lambda)) {
    # Model was estimated using multiple lambdas
    # but user only wants to use one of them
    Phi_hat <- Phi_hat[, , which(mod$lambdas == lambda)]
    phi_0 <- mod$phi0hat[, , which(mod$lambdas == lambda)]
  }

  fcst <- NULL
  if (!is.matrix(Phi_hat) & is.null(lambda)){
    # Forecasting for all lambdas
    fcst <- array(dim = c(h, mod$k, dim(Phi_hat)[3]))
    for (slice in 1:dim(Phi_hat)[3]){
      fcst[, , slice] <- .recursiveforecast(mod$Y, Phi_hat[, , slice],
                                            phi_0[, , slice], init, h)
    }
  }
  else {
    fcst <- .recursiveforecast(mod$Y, Phi_hat, phi_0, init, h)
  }

  # Checking whether cv was performed. In this case mod$lambdas > 1 but only
  # one was truly used
  if (!is.na(mod$lambda_SEopt)) lambda <- mod$lambda_SEopt
  if (is.null(lambda)) lambda <- mod$lambdas

  out <- list(
    fcst = fcst,
    h = h,
    lambda = lambda,
    Y = mod$Y
  )

  class(out) <- "bigtime.recursiveforecast"
  out
}

#' Recursively forecast a VAR
#'
#' This function is not meant to be directly called by the user
#' and is only a helper function
#'
#' @param Y data matrix
#' @param Phi_hat coefficient matrix
#' @param phi_0 constent terms
#' @param init starting values of recursive forecast. Usually the last p
#' observations of each variable. Must be in the form
#' (y_1t, y_2t, y_1t-1, y_2t-1, ...)
#' @param h forecast horizon
#' @return Returns a matrix of forecasts
.recursiveforecast <- function(Y, Phi_hat, phi_0, init, h){
  # Constructing the companion form for quicker forecasting
  k <- nrow(Phi_hat)
  p <- ncol(Phi_hat)/k

  I <- diag(1, k*(p-1), k*(p-1))
  O <- matrix(0, k*(p-1), k)

  FF <- rbind(Phi_hat, cbind(I, O))
  const <- c(phi_0, rep(0, k*(p-1)))

  # Getting the starting values, e.g. the last p observations of each series
  yy <- apply(Y, 2, rev)
  init <- c(t(yy))[1:(k*p)]

  # Forecasting
  fcst <- recursiveforecast_cpp(init, FF, const, h)
  fcst <- fcst[-1, 1:k, drop = FALSE]
  colnames(fcst) <- colnames(Y)
  fcst
}


#' Plots Recursive Forecasts
#'
#' Plots the recursive forecast obtained using bigtime::recursiveforecast
#' When forecasts were made for multiple lambdas and lmbda is not a single
#' number, then a ribbon will be plotted that reaches from the minimum estimate
#' of all lambdas to the maximum.
#'
#' If lmbda is of length one or forecsts were made using only one lambda,
#' then only a line will be plotted.
#'
#' Default names for series are Y1, Y2, ... if the original data does
#' not have any column names
#'
#' @param fcst Recursive Forecast obtained using bigtime::recursiveforecast
#' @param series Series name. If original data has not names, then use Y1 for
#' the first series, Y2 for the second, and so on
#' @param lmbda lambdas to be used for plotting. If forecast was done using only
#' one lambda, then this will be ignored.
#' @param last_n last n observations of the original data to include in the plot
#' @export
#' @return Returns a ggplot
plot.bigtime.recursiveforecast <- function(fcst, series, lmbda=NULL,
                                           last_n = floor(nrow(fcst$Y)*0.1)){
  # Setting up Y so it is in the right form
  Y <- as.data.frame(fcst$Y)
  cnames <- colnames(Y)
  if (is.null(cnames)) cnames <- paste0("Y", 1:dim(Y)[2])
  colnames(Y) <- cnames

  # Checking whether the forecast was done using multiple lambdas
  if (length(fcst$lambda) > 1){
    # We forcasted for multiple lambdas
    fcast <- reduce_cube(fcst$fcst, fcst$lambda, "lambda")
  }else {
    fcast <- as.data.frame(fcst$fcst)
    if (is.null(colnames(fcst$fcst))) colnames(fcast) <- paste0("Y", 1:dim(fcast)[2])
    if (is.null(colnames(fcast))) colnames(fcast) <- paste0("Y", 1:dim(fcast)[2])
    lambda <- rep(fcst$lambda, dim(fcast)[1])
    fcast <- cbind(fcast, lambda)
  }

  if(!require(tidyverse, quietly = TRUE)) stop("tidyverse must be installed")
  Y <- as_tibble(Y[, series, drop = FALSE]) %>% mutate(t = 1:n())
  fcast <- as_tibble(fcast[, c(series, "lambda")]) %>%
    group_by(lambda) %>%
    mutate(t = 1:n()+nrow(Y)) %>%
    ungroup()

  if (length(fcst$lambda) == 1 | length(lmbda) == 1){
    # The model was either forecasted using only one lambda
    # or we wish to plot only for one lambda
    if (is.null(lmbda)) lmbda <- fcst$lambda
    if (!(lmbda %in% fcast$lambda)) stop("Selected lambda was not used for forecasting.")
    fcast <- fcast %>% filter(lambda == lmbda)

    p <- Y %>%
      filter(t >= n() - last_n) %>%
      ggplot() +
      geom_line(aes_string("t", series), color = "black") +
      geom_line(data = fcast, aes_string("t", series, group = "lambda"), color = "coral3") +
      theme_bw()
  }else {
    # Plotting for multiple lambdas
    # In this case we plot a ribbon
    message("Plotting ribbon from min to max because multiple lambdas are used")
    if (!is.null(lmbda)) fcast <- fcast %>% filter(lambda %in% lmbda)
    fcast <- fcast %>%
      group_by(t) %>%
      mutate_at(.vars = vars(!!series), .funs = list(min = min, max = max)) %>%
      summarise_all(.funs = function(x) x[[length(x)]])

    p <- Y %>%
      filter(t >= n() - last_n) %>%
      ggplot() +
      geom_line(aes_string("t", series), color = "black") +
      geom_ribbon(data = fcast, aes(x = t, ymin = min, ymax = max), alpha = 0.5, fill = "coral3", color = "coral3") +
      theme_bw()
  }

  p + xlab("Time")
}

#' Checks whether a VAR is stable
#'
#' Using a model estimated by bigtime::sparseVAR, this function checks whether
#' the resulting VAR is stable
#'
#' @param mod model estimated using bigtime::sparseVAR. Can only be a model
#' with one coefficient vector. Hence, the model must either be estimated using
#' cv=TRUE or by giving a single lambda value
#' @param verbose If TRUE, then the actual maximum eigenvalue of the companion
#' matrix will be printed to the console. Default is FALSE
#' @export
#' @return Returns TRUE if the VAR is stable and FALSE otherwise
is.stable <- function(mod, verbose = FALSE){
  if (!("bigtime.VAR") %in% class(mod)) stop("Model is not a VAR model estimated using bigtime::sparseVAR")
  if (!is.matrix(mod$Phihat)) stop("Model contains multiple coefficient estimates. It is not clear which model is meant. Please reduce the model to one coefficient estimate")
  Phi_hat <- mod$Phihat
  k <- mod$k
  p <- mod$p
  I <- diag(1, k*(p-1), k*(p-1))
  O <- matrix(0, k*(p-1), k)

  FF <- rbind(Phi_hat, cbind(I, O))
  max_eigval <- max(abs(eigen(FF)$value))
  if (verbose) cat("Maximum eigenvalu of Companion Matrix: ", max_eigval, "\n")
  if (max_eigval < 1) return(TRUE)
  FALSE
}


#' Reduces the third dimension of a cube and returns a data frame
#'
#' This is only meant to be a helper function and not meant to be called by
#' the user itself
#'
#' @param cube some 3D array
#' @param dim3_names vector of length dim(cube)[3] containing a name for
#' each slice
#' @param name what should the column containing the dim3_names be called?
#' @return Returns a data frame contructed from the cube
reduce_cube <- function(cube, dim3_names, name = deparse(substitute(dim3_names))){
  cnames <- colnames(cube)
  if (is.null(cnames)) {
    message("Cube has no column names. Defaulting to Y...")
    cnames <- paste0("Y", 1:dim(cube)[2])
  }
  out <- NULL
  for (slice in 1:dim(cube)[3]){
    df <- data.frame(cube[, , slice])
    colnames(df) <- cnames
    des <- rep(dim3_names[slice], dim(cube)[1])
    df <- cbind(df, des)
    colnames(df)[[dim(cube)[2] + 1]] <- name
    if (is.null(out)) out <- df
    else out <- rbind(out, df)
  }
  out
}









