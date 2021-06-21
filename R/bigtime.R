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
#' The function \link{recursiveforecast} can be used to recursively forecast a VAR model.
#' The function \link{is.stable} returns whether an estimated VAR model is stable.
#' The function \link{diagnostics_plot} returns a plot of the fitted vs. observed values as well as of the residuals. 
#' The functions \link{fitted} and \link{residuals} return the fitted, respectively the residuals of the estimated model.
#' The function \link{simVAR} can be used to simulate a VAR model with various sparsity patterns.
#'
#' @author Ines Wilms <i.wilms@maastrichtuniversity.nl>, Jacob Bien, David S. Matteson, Sumanta Basu, Will Nicholson, Enrico Wegner
#' @references Nicholson William B., Wilms Ines, Bien Jacob and Matteson David S. (2020), “High-dimensional forecasting via
#' interpretable vector autoregression”, Journal of Machine Learning Research, 21(166), 1-52.
#' @references Wilms Ines, Sumanta Basu, Bien Jacob and Matteson David S. (2021), “Sparse Identification and
#' Estimation of Large-Scale Vector AutoRegressive Moving Averages”, Journal of the American Statistical Association, doi: 10.1080/01621459.2021.1942013.

#' @docType package
#' @name bigtime
#' @examples
#' # Fit a sparse VAR model
#' data(var.example)
#' VARfit <- sparseVAR(scale(Y.var), selection = "cv") # sparse VAR using time series cross-validation
#' Lhat <- lagmatrix(fit=VARfit) # get estimated lagmatrix
#' VARforecast <- directforecast(fit=VARfit, h=1) # get one-step ahead forecasts
NULL

