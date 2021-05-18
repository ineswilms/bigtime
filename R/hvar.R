#' Sparse Estimation of the Vector AutoRegressive (VAR) Model
#' @param Y A \eqn{T} by \eqn{k} matrix of time series. If k=1, a univariate autoregressive model is estimated.
#' @param p User-specified maximum autoregressive lag order of the VAR. Typical usage is to have the program compute its own maximum lag order based on the time series length.
#' @param h Desired forecast horizon in time-series cross-validation procedure.
#' @param VARlseq User-specified grid of values for regularization parameter corresponding to sparse penalty. Typical usage is to have the program compute
#' its own grid. Supplying a grid of values overrides this. WARNING: use with care.
#' @param VARgran User-specified vector of granularity specifications for the penalty parameter grid:  First element specifies
#' how deep the grid should be constructed. Second element specifies how many values the grid should contain.
#' @param cvcut Proportion of observations used for model estimation in the time series cross-validation procedure. The remainder is used for forecast evaluation.
#' @param eps a small positive numeric value giving the tolerance for convergence in the proximal gradient algorithm.
#' @param VARpen "HLag" (hierarchical sparse penalty) or "L1" (standard lasso penalty) penalization.
#' @param cv Logical, whether time-series cross-validation needs to be performed (TRUE) or not (FALSE) for selecting the sparsity parameter. If cv=FALSE, the argument cvcut is redundant.
#' @export
#' @return A list with the following components
#' \item{Y}{\eqn{T} by \eqn{k} matrix of time series.}
#' \item{k}{Number of time series.}
#' \item{p}{Maximum autoregressive lag order of the VAR.}
#' \item{Phihat}{Matrix of estimated autoregressive coefficients of the VAR.}
#' \item{phi0hat}{vector of VAR intercepts.}
#' \item{series_names}{names of time series}
#' \item{lambdas}{sparsity parameter grid}
#' \item{MSFEcv}{MSFE cross-validation scores for each value of the sparsity parameter in the considered grid}
#' \item{MSFEcv_all}{MSFE cross-validation full output}
#' \item{lambda_opt}{Optimal value of the sparsity parameter as selected by the time-series cross-validation procedure}
#' \item{lambda_SEopt}{Optimal value of the sparsity parameter as selected by the time-series cross-validation
#'       procedure and after applying the one-standard-error rule.  This is the value used.}
#' \item{h}{Forecast horizon h}
#' @references Nicholson William B., Bien Jacob and Matteson David S. (2017), "High-dimensional forecasting via interpretable vector autoregression," Journal of Machine Learning Research, 21(166), 1-52.
#' @seealso \link{lagmatrix} and \link{directforecast}
#' @examples
#' data(Y)
#' VARfit <- sparseVAR(Y) # sparse VAR
#' ARfit <- sparseVAR(Y[,2]) # sparse AR
sparseVAR <- function(Y, p=NULL, VARpen="HLag", VARlseq=NULL, VARgran=NULL,
                      cv = FALSE, cvcut=0.9, h=1,  eps=1e-3){

  # Check Inputs
  if(!is.matrix(Y)){

    if(is.vector(Y) & length(Y)>1){
      Y <- matrix(Y, ncol=1)
    }else{
      stop("Y needs to be a matrix of dimension T by k")
    }

  }


  if(nrow(Y)<10){
    stop("The time series length is too small.")
  }

  if(!is.null(colnames(Y))){ # time series names
    series_names <- colnames(Y)
  } else {
    series_names <- NULL
  }

  if(h<=0){
    stop("The forecast horizon h needs to be a strictly positive integer")
  }

  if(!is.null(p)){
    if(p<=0){
      stop("The maximum autoregressive order of VAR needs to be a strictly positive integer")
    }
  }

  if((!is.vector(VARlseq) & !is.null(VARlseq)) ){ # | length(VARlseq)==1
    stop("The regularization parameter VARlseq needs to be a vector of length=>1 or NULL otherwise")
  }

  if(any((VARgran<=0)==T)){
    stop("The granularity parameters needs to be a strictly positive integers")
  }

  if(!((cvcut<1) & (cvcut>0))){
    stop("cvcut needs to be a number between 0 and 1")
  }


  if(!is.element(VARpen, c("HLag", "L1"))){
    stop("The type of penalization VARpen needs to be either HLag (hierarchical sparse penalization) or L1 (standard lasso penalization)")
  }

  if(!is.null(p)){
    if(VARpen=="HLag" & p<=1){
      stop("HLag penalization is only supported for p>1. Use L1 as VARpen instead")
    }
  }


  if(eps<=0){
    stop("The convergence tolerance parameter eps needs to be a small positive number")
  }


  if(is.null(p)){
    p <- floor(1.5*sqrt(nrow(Y)))
  }

  # Regularization grid
  if(is.null(VARgran)){
    VARgran1 <- 10^2
    VARgran2 <- 10
  }else{
    VARgran1 <- VARgran[1]
    VARgran2 <- VARgran[2]
  }

  if(!is.null(VARlseq)){
    VARgran1 <- max(VARlseq)/min(VARlseq)
    VARgran2 <- length(VARlseq)
  }

  if(length(VARlseq)==1){ # If user only provides one lambda value, don't do cross-validation
    cv <- FALSE
    warning("No cross-validation is performed since only one sparsity parameter is provided.")
  }

  # Set the grid of sparsity parameters
  VARdata <- HVARmodel(Y=Y, p=p, h=h)
  k <- VARdata$k # Number of time series

  if(cv){ # Time series cross-validation to get optimal sparsity parameter
    VARcv <- HVAR_cv(Y=Y, p=p, h=h, lambdaPhiseq=VARlseq, gran1=VARgran1, gran2=VARgran2, T1.cutoff=cvcut, eps=eps, type=VARpen)

    # Var estimation with selected regularization parameter
    VARmodel <- HVAR(fullY=VARdata$fullY, fullZ=VARdata$fullZ, p=VARdata$p, k=VARdata$k, lambdaPhi=VARcv$lambda_opt_oneSE, eps=eps, type=VARpen)

  }else{ # No time series cross-validation

    Phis <- array(NA, c(k, k*p, VARgran2)) # Estimates AR coefficients for each value in the grid
    phi0s <- array(NA, c(k, 1, VARgran2)) # Estimates of constants for each value in the grid

    fullY <- VARdata$fullY # response matrix

    if(k==1){
      fullY <- matrix(fullY, ncol=1)
    }
    fullZ <- VARdata$fullZ # design matrix


    # Get lambda grid if not specified by the user
    if(is.null(VARlseq)){
      jj <- .lfunction3(p, k)

      if(VARpen=="HLag"){
        VARlseq <- .LambdaGridE(VARgran1, VARgran2, jj, fullY, fullZ, "HVARELEM", p, k,
                                MN=F, alpha=1/(k+1), C=rep(1,p))
      }

      if(VARpen=="L1"){
        VARlseq <- .LambdaGridE(VARgran1, VARgran2, jj, fullY, fullZ,"Basic",p,k,MN=F,alpha=1/(k+1),C=rep(1,p))
      }

    }

    for(il in 1:length(VARlseq)){ # For now in R
      VARmodel <- HVAR(fullY=VARdata$fullY, fullZ=VARdata$fullZ, p=VARdata$p, k=VARdata$k, lambdaPhi=VARlseq[il], eps=eps, type=VARpen)
      Phis[,,il] <- VARmodel$Phi
      phi0s[,,il] <- VARmodel$phi
    }

  }

  if(cv){
    out <- list("k"=k, "Y"=Y, "p"=p, "Phihat"=VARmodel$Phi, "phi0hat"=VARmodel$phi,
                "series_names"=series_names, "lambdas"=VARcv$lambda,
                "MSFEcv"=VARcv$MSFE_avg, "MSFE_all"=VARcv$MSFE_all,
                "lambda_SEopt"=VARcv$lambda_opt_oneSE,"lambda_opt"=VARcv$lambda_opt, "h"=h)
  }else{
    out <- list("k"=k, "Y"=Y, "p"=p, "Phihat"=Phis, "phi0hat"=phi0s,
                "series_names"=series_names, "lambdas"=VARlseq,
                "MSFEcv"=NA, "MSFE_all"=NA,
                "lambda_SEopt"=NA,"lambda_opt"=NA, "h"=h)
  }

  class(out) <- "bigtime.VAR"
  out
}

HVARmodel<-function(Y, p, h=1){
  # Preliminaries
  k <- ncol(Y) # Number of Endogenous variables

  # Lagged predictor matrices
  DATAY <- embed(Y, dimension=p+h) # collect data to compute h-step ahead forecasts when selecting lambda
  fullY <- DATAY[, 1:k] # in its columns Y1, Y2, ... Yq
  fullZ <- t(DATAY[, (ncol(DATAY)-k*p+1):ncol(DATAY)])

  out<-list("fullY"=fullY, "fullZ"=fullZ, "k"=k, "p"=p)
}

HVAR<-function(fullY, fullZ, p, k, lambdaPhi, eps=1e-5, type="HLag"){

  # HVAR estimates
  if(type=="HLag"){
    if(k==1){
      fullY <- matrix(fullY, ncol=1)
    }

    Phi <- HVARElemAlgcpp(array(0, dim=c(k,k*p+1,1)), fullY, fullZ, lambdaPhi, eps, p)
  }

  # L1 estimates
  if(type=="L1"){
    Phi <- lassoVARFistcpp(array(0,dim=c(k,k*p+1,1)), fullY, fullZ, lambdaPhi, eps, p)
  }

  resids <- t(t(fullY) - Phi[,,1]%*%rbind(rep(1, ncol(fullZ)), fullZ))
  Yhat <- t(Phi[,, 1]%*%rbind(rep(1, ncol(fullZ)), fullZ))

  # Output
  out <- list("Phi"=Phi[,-1,1], "phi0"=Phi[,1,1], "resids"=resids, "Yhat"=Yhat, "Y"=fullY)

}

HVAR_cv<-function(Y, p, h=1, lambdaPhiseq=NULL, gran1 = 10^2, gran2=10, T1.cutoff=0.9, eps=1e-5, type="HLag"){

  # Get response and predictor matrix
  HVARmodelFIT <- HVARmodel(Y=Y, p=p, h=h)
  k <- HVARmodelFIT$k # Number of time series
  fullY <- HVARmodelFIT$fullY # response matrix
  if(k==1){
    fullY <- matrix(fullY, ncol=1)
  }
  fullZ <- HVARmodelFIT$fullZ # design matrix


  # Get lambda grid if not specified by the user
  if(is.null(lambdaPhiseq)){
    jj <- .lfunction3(p, k)

    if(type=="HLag"){
      lambdaPhiseq <- .LambdaGridE(gran1, gran2, jj, fullY, fullZ,"HVARELEM", p, k,
                                       MN=F, alpha=1/(k+1), C=rep(1,p))
    }

    if(type=="L1"){
      lambdaPhiseq <- .LambdaGridE(gran1, gran2, jj, fullY, fullZ,"Basic",p,k,MN=F,alpha=1/(k+1),C=rep(1,p))
    }

  }

  # Time Series cross-validation loop
  # Choose sparsity parameter that minimizes h-step ahead mean squared prediction error.
  n <- nrow(fullY)
  T1 <- floor(T1.cutoff*n)
  tseq <- T1:(n-1)

  estim <- (type=="HLag")*2 + (type=="L1")*1
  # Time-series cross-validation for tuning parameter selection
  my_cv <- HVAR_cvaux_loop_cpp(Y = fullY, Z = fullZ, tseq = tseq, gamm = lambdaPhiseq,  eps = eps, p = p, estim = estim)
  MSFEmatrix <- my_cv$MSFEcv
  rownames(MSFEmatrix) <- paste0("t=", tseq)
  colnames(MSFEmatrix) <- paste0("lambda=", lambdaPhiseq)
  sparsitymatrix <- my_cv$sparsitycv
  colnames(sparsitymatrix) <- colnames(MSFEmatrix)
  rownames(sparsitymatrix) <- rownames(MSFEmatrix)
  MSFE_avg <- apply(MSFEmatrix, 2, mean)
  lambda_opt <- lambdaPhiseq[which.min(MSFE_avg)]

  if(length(lambda_opt)==0){
    lambda_opt <- lambdaPhiseq[round(median(1:length(lambdaPhiseq)))]
  }else{
    if(is.na(lambda_opt)){
      lambda_opt <- lambdaPhiseq[round(median(1:length(lambdaPhiseq)))]
    }
  }

  # One-standard error rule to determine optimal lambda values
  lambda_opt_oneSE <- NA
  lambda_optinit <- lambda_opt
  MSFE_sd <- apply(MSFEmatrix, 2, sd)/sqrt(nrow(MSFEmatrix))
  MSFE_flagOK <- MSFE_avg < min(MSFE_avg, na.rm=T) + MSFE_sd[which.min(MSFE_avg)]
  MSFEnew <- MSFE_avg
  MSFEnew[!MSFE_flagOK] <- NA
  sparsitynew <- apply(sparsitymatrix,2,mean)
  sparsitynew[!MSFE_flagOK] <- NA

  lambda_opt_oneSE <- lambdaPhiseq[which.min(sparsitynew)]

  if(length(lambda_opt_oneSE)==0){
    lambda_opt_oneSE <- lambda_optinit
  }else{
    if(is.na(lambda_opt_oneSE)){
      lambda_opt_oneSE <- lambda_optinit
    }
  }

  gridflag_oneSE <- lambda_opt_oneSE==min(lambdaPhiseq)
  if(lambda_opt_oneSE==min(lambdaPhiseq)){
    warning("Lower bound of lambda grid is selected: higher value of first granularity parameter is recommended as an input")
  }


  # Output
  out <- list("lambda" = lambdaPhiseq, "lambda_opt_oneSE" = lambda_opt_oneSE, "MSFE_all" = MSFEmatrix, "MSFE_avg" = MSFE_avg, "flag_oneSE" = gridflag_oneSE,
              "lambda_opt" = lambda_opt)

}
