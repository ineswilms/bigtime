#' bigtime: A package for obtaining sparse estimates of large time series models.
#'
#' The bigtime package provides sparse estimators for three large time series models:
#' Vector AutoRegressive Models,
#' Vector AutoRegressive Models with Exogenous variables,
#' and Vector AutoRegressive Moving Average Models. The univariate cases are also supported.
#'
#' @details To use the facilities of this package, start with a T by k time series matrix Y (for the VAR and VARMA), and an exogenous time series matrix X (for the VARX).
#' Run \link{sparseVAR}, \link{sparseVARX} or \link{sparseVARMA} to get the estimated model.
#' The function \link{lagmatrix} returns the lag matrix of estimated coefficients of the estimated model.
#' The function \link{directforecast} gives h-step ahead forecasts based on the estimated model.
#'
#' @author Ines Wilms <ines.wilms@kuleuven.be>, Jacob Bien, David Matteson, Sumanta Basu
#' @references Nicholson William B., Bien Jacob and Matteson David S. (2017), "High Dimensional Forecasting via Interpretable Vector Autoregression"
#' arXiv preprint arXiv:1412.5250v2.
#' @references Wilms Ines, Sumanta Basu, Bien Jacob and Matteson David S. (2017), "Sparse Identification and Estimation of High-Dimensional Vector AutoRegressive Moving Averages"
#' arXiv preprint arXiv:1707.09208.
#' @references Wilms Ines, Sumanta Basu, Bien Jacob and Matteson David S. (2017), "Interpretable Vector AutoRegressions with
#' Exogenous Time Series" arXiv preprint.
#'
#' @docType package
#' @name bigtime
#' @examples
#' # Fit a sparse VAR model
#' data(Y)
#' VARfit <- sparseVAR(Y) # sparse VAR
#' Lhat <- lagmatrix(fit=VARfit, model="var") # get estimated lagmatrix
#' VARforecast <- directforecast(fit=VARfit, model="VAR", h=1) # get one-step ahead forecasts
NULL
