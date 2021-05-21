


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
#' @example
#' sim_data <- bigtime::simVAR(200, 5, 5)
#' summary(sim_data)
#' mod <- bigtime::sparseVAR(sim_data$Y)
#' bigtime::is.stable(mod)
#' fcst_recursive <- bigtime::recursiveforecast(mod, h = 4)
#' plot(fcst_recursive, series = "Y1")
#' fcst_direct <- bigtime::directforecast(mod, "VAR")
#' fcst_direct
#' fcst_recursive$fcst
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
                                            phi_0[, , slice], h)
    }
  }
  else {
    fcst <- .recursiveforecast(mod$Y, Phi_hat, phi_0, h)
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
#' @param h forecast horizon
#' @return Returns a matrix of forecasts
.recursiveforecast <- function(Y, Phi_hat, phi_0, h){
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
#' @param x Recursive Forecast obtained using bigtime::recursiveforecast
#' @param series Series name. If original data has not names, then use Y1 for
#' the first series, Y2 for the second, and so on
#' @param lmbda lambdas to be used for plotting. If forecast was done using only
#' one lambda, then this will be ignored.
#' @param last_n last n observations of the original data to include in the plot
#' @param ... Not currently used
#' @export
#' @return Returns a ggplot
plot.bigtime.recursiveforecast <- function(x, series=NULL, lmbda=NULL,
                                           last_n = floor(nrow(fcst$Y)*0.1), ...){
  cn <- colnames(x$fcst)
  if (is.na(series)) series <- ifelse(is.null(cn), "Y1", cn[[1]])
  fcst <- x
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

  # if(!require(tidyverse, quietly = TRUE)) stop("tidyverse must be installed")
  Y <- dplyr::as_tibble(Y[, series, drop = FALSE]) %>% dplyr::mutate(t = 1:dplyr::n())
  fcast <- dplyr::as_tibble(fcast[, c(series, "lambda")]) %>%
    dplyr::group_by(lambda) %>%
    dplyr::mutate(t = 1:dplyr::n()+nrow(Y)) %>%
    dplyr::ungroup()

  if (length(fcst$lambda) == 1 | length(lmbda) == 1){
    # The model was either forecasted using only one lambda
    # or we wish to plot only for one lambda
    if (is.null(lmbda)) lmbda <- fcst$lambda
    if (!(lmbda %in% fcast$lambda)) stop("Selected lambda was not used for forecasting.")
    fcast <- fcast %>% dplyr::filter(lambda == lmbda)

    p <- Y %>%
      dplyr::filter(t >= dplyr::n() - last_n) %>%
      ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes_string("t", series), color = "black") +
      ggplot2::geom_line(data = fcast, ggplot2::aes_string("t", series, group = "lambda"), color = "coral3") +
      ggplot2::theme_bw()
  }else {
    # Plotting for multiple lambdas
    # In this case we plot a ribbon
    message("Plotting ribbon from min to max because multiple lambdas are used")
    if (!is.null(lmbda)) fcast <- fcast %>% dplyr::filter(lambda %in% lmbda)
    fcast <- fcast %>%
      dplyr::group_by(t) %>%
      dplyr::mutate_at(.vars = dplyr::vars(!!series), .funs = list(min = min, max = max)) %>%
      dplyr::summarise_all(.funs = function(x) x[[length(x)]])

    p <- Y %>%
      dplyr::filter(t >= dplyr::n() - last_n) %>%
      ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes_string("t", series), color = "black") +
      ggplot2::geom_ribbon(data = fcast, ggplot2::aes(x = t, ymin = min, ymax = max), alpha = 0.5, fill = "coral3", color = "coral3") +
      ggplot2::theme_bw()
  }

  p + ggplot2::xlab("Time")
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









