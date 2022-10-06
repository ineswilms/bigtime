#' Sparse Estimation of the Vector AutoRegressive  with Exogenous Variables X (VARX) Model
#' @param Y A \eqn{T} by \eqn{k} matrix of time series. If k=1, a univariate autoregressive model is estimated.
#' @param X A \eqn{T} by \eqn{m} matrix of time series.
#' @param p User-specified maximum endogenous autoregressive lag order. Typical usage is to have the program compute its own maximum lag order based on the time series length.
#' @param s User-specified maximum exogenous autoregressive lag order. Typical usage is to have the program compute its own maximum lag order based on the time series length.
#' @param h Desired forecast horizon in time-series cross-validation procedure.
#' @param VARXlPhiseq User-specified grid of values for regularization parameter corresponding to the endogenous autoregressive coefficients in the VARX. Typical usage is to have the program compute
#' its own grid. Supplying a grid of values overrides this. WARNING: use with care.
#' @param VARXPhigran User-specified vector of granularity specifications for the penalty parameter grid corresponding to the endogenous autoregressive coefficients in the VARX:  First element specifies
#' how deep the grid should be constructed. Second element specifies how many values the grid should contain.
#' @param VARXlBseq User-specified grid of values for regularization parameter corresponding to the exogenous autoregressive coefficients in the VARX. Typical usage is to have the program compute
#' its own grid. Supplying a grid of values overrides this. WARNING: use with care.
#' @param VARXBgran User-specified vector of granularity specifications for the penalty parameter grid corresponding to the exogenous autoregressive coefficients in the VARX:  First element specifies
#' how deep the grid should be constructed. Second element specifies how many values the grid should contain.
#' @param eps a small positive numeric value giving the tolerance for convergence in the proximal gradient algorithm.
#' @param cvcut Proportion of observations used for model estimation in the time series cross-validation procedure. The remainder is used for forecast evaluation.
#' @param VARXalpha a small positive regularization parameter value corresponding to squared Frobenius penalty. The default is zero.
#' @param VARXpen "HLag" (hierarchical sparse penalty) or "L1" (standard lasso penalty) penalization in VARX.
#' @param selection Model selection method to be used. Default is none, which will return all values for all penalisations.
#' @param check_std Check whether data is standardised. Default is TRUE and is not recommended to be changed
#' @param verbose Logical to print value of information criteria for each lambda together with selection. Default is FALSE
#' @export
#' @return A list with the following components
#' \item{Y}{\eqn{T} by \eqn{k} matrix of endogenous time series.}
#' \item{X}{\eqn{T} by \eqn{m} matrix of exogenous time series.}
#' \item{k}{Number of endogenous time series.}
#' \item{m}{Number of exogenous time series.}
#' \item{p}{Maximum endogenous autoregressive lag order of the VARX.}
#' \item{s}{Maximum exogenouss autoregressive lag order of the VARX.}
#' \item{Phihat}{Matrix of estimated endogenous autoregressive coefficients.}
#' \item{Bhat}{Matrix of estimated exogenous autoregressive coefficients.}
#' \item{phi0hat}{vector of VARX intercepts.}
#' \item{exogenous_series_names}{names of the exogenous time series}
#' \item{endogenous_series_names}{names of the endogenous time series}
#' \item{lambdaPhi}{sparsity parameter grid corresponding to endogenous autoregressive parameters}
#' \item{lambdaB}{sparsity parameter grid corresponding to exogenous autoregressive parameters}
#' \item{lambdaPhi_opt}{Optimal value of the sparsity parameter (corresponding to the endogenous autoregressive parameters) as selected by the time-series cross-validation procedure}
#' \item{lambdaPhi_SEopt}{Optimal value of the sparsity parameter (corresponding to the endogenous autoregressive parameters) as selected by the time-series cross-validation procedure and after applying the one-standard-error rule}
#' \item{lambdaB_opt}{Optimal value of the sparsity parameter (corresponding to the exogenous autoregressive parameters) as selected by the time-series cross-validation procedure}
#' \item{lambdaB_SEopt}{Optimal value of the sparsity parameter (corresponding to the exogenous autoregressive parameters) as selected by the time-series cross-validation procedure and after applying the one-standard-error rule}
#' \item{MSFEcv}{MSFE cross-validation scores for each value in the two-dimensional sparsity grid}
#' \item{h}{Forecast horizon h}
#' @references Wilms Ines, Sumanta Basu, Bien Jacob and Matteson David S. (2017), “Interpretable vector autoregressions with exogenous time series”, NIPS 2017 Symposium on Interpretable Machine Learning, arXiv:1711.03623.
#' @seealso \link{lagmatrix} and \link{directforecast}
#' @examples
#' data(varx.example)
#' VARXfit <- sparseVARX(Y=scale(Y.varx), X=scale(X.varx)) # sparse VARX
#' y <- matrix(Y.varx[,1], ncol=1)
#' ARXfit <- sparseVARX(Y=y, X=X.varx) # sparse ARX
sparseVARX <- function(Y, X, p=NULL, s=NULL, VARXpen="HLag", VARXlPhiseq=NULL, VARXPhigran=NULL,
                       VARXlBseq=NULL,  VARXBgran=NULL, VARXalpha=0, h=1, cvcut=0.9, eps=10^-3,
                       selection = c("none", "cv", "bic", "aic", "hq"),
                       check_std = TRUE, verbose = FALSE){

  ######################################
  #### Check Inputs and Preparation ####
  ######################################
  selection <- match.arg(selection)
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
  if(!is.matrix(X)){
    if(is.vector(X) & length(X)>1){
      X <- matrix(X, ncol=1)
    }else{
      stop("X needs to be a matrix of dimension T by m")
    }
  }
  if(nrow(X)!=nrow(Y)){
    stop("Y and X need to have the same number of rows (observations).")
  }
  if(!is.null(colnames(Y))){ # time series names
    endogenous_series_names <- colnames(Y)
  } else {
    endogenous_series_names <- NULL
  }
  if(!is.null(colnames(X))){ # time series names
    exogenous_series_names <- colnames(X)
  } else {
    exogenous_series_names <- NULL
  }
  if(!is.null(p)){
    if(p<=0){
      stop("The maximum endogenous autoregressive order needs to be a strictly positive integer")
    }
  }
  if(!is.null(s)){
    if(s<=0){
      stop("The maximum exogenous autoregressive order needs to be a strictly positive integer")
    }
  }
  if(h<=0){
    stop("The forecast horizon h needs to be a strictly positive integer")
  }
  if(!((cvcut<1) & (cvcut>0))){
    stop("cvcut needs to be a number between 0 and 1")
  }
  if( (!is.vector(VARXlPhiseq) & !is.null(VARXlPhiseq))){
    stop("The regularization parameter grid VARXlPhiseq needs to be a vector of length >= 1 or NULL otherwise")
  }
  if (!is.null(VARXlPhiseq) & length(VARXlPhiseq)==1 & selection!="none"){
    if (!is.null(VARXlBseq) & length(VARXlBseq)==1){
      stop("When providing a single penalisation parameter, selection must be 'none'")
    }
  }
  if(any((VARXPhigran<=0)==T)){
    stop("The granularity parameters need to be a strictly positive integer")
  }
  if((!is.vector(VARXlBseq) & !is.null(VARXlBseq))){
    stop("The regularization parameter VARXlBseq needs to be a vector of length >=1 or NULL otherwise")
  }
  if (!is.null(VARXlBseq) & length(VARXlBseq)==1 & selection!="none"){
    if (!is.null(VARXlPhiseq) & length(VARXlPhiseq)==1){
      stop("When providing a single penalisation parameter, selection must be 'none'")
    }
  }
  if(any((VARXBgran<=0)==T)){
    stop("The granularity parameters need to be a strictly positive integer")
  }
  if(VARXalpha<0){
    stop("The regularization paramter VARXalpha needs to be a equal to zero or a small positive number")
  }
  if(!is.element(VARXpen, c("HLag", "L1"))){
    stop("The type of penalization VARXpen needs to be either HLag (hierarchical sparse penalization) or L1 (standard lasso penalization)")
  }
  if(!is.null(p) | !is.null(s)){
    if(is.null(p)){
      p <- floor(1.5*sqrt(nrow(Y)))
    }
    if(is.null(s)){
      s <- floor(1.5*sqrt(nrow(Y)))
    }
    if(VARXpen=="HLag" & (p<=1 |s<=1)){
      stop("HLag penalization in VARX is only supported for p and s larger than 1. Use L1 as VARXpen instead")
    }
  }
  if(eps<=0){
    stop("The convergence tolerance parameter eps needs to be a small positive number")
  }
  if (check_std) .check_if_standardised(Y)
  if (check_std) .check_if_standardised(X)

  # Set maximum orders
  if(is.null(p)){
    p <- floor(1.5*sqrt(nrow(Y)))
  }
  if(is.null(s)){
    s <- floor(1.5*sqrt(nrow(Y)))
  }

  # Regularization grid
  if(is.null(VARXPhigran)){
    VARXPhigran1 <- 10^2
    VARXPhigran2 <- 10
  }else{
    VARXPhigran1 <- VARXPhigran[1]
    VARXPhigran2 <- VARXPhigran[2]
  }
  if(!is.null(VARXlPhiseq)){
    VARXPhigran1 <- max(VARXlPhiseq)/min(VARXlPhiseq)
    VARXPhigran2 <- length(VARXlPhiseq)
  }
  if(is.null(VARXBgran)){
    VARXBgran1 <- 10^2
    VARXBgran2 <- 10
  }else{
    VARXBgran1 <- VARXBgran[1]
    VARXBgran2 <- VARXBgran[2]
  }
  if(!is.null(VARXlBseq)){
    VARXBgran1 <- max(VARXlBseq)/min(VARXlBseq)
    VARXBgran2 <- length(VARXlBseq)
  }

  # Preparing VARX data
  VARXdata <- HVARXmodel(Y=Y, X=X, p=p, s=s, h=h)
  estim <- (VARXpen=="HLag")*2 + (VARXpen=="L1")*1
  if(VARXdata$k==1){
    VARXdata$fullY <- as.matrix(VARXdata$fullY)
  }
  k <- ncol(Y)
  m <- ncol(X)
  Phihat <- NULL
  phi0hat <- NULL
  Bhat <- NULL

  ############################
  #### END Check and Prep ####
  ############################

  if (selection == "cv"){
    # Get optimal sparsity parameter via time series cross-validation
    VARXcv <- HVARX_cv(Y=Y, X=X, p=p, s=s, h=h, lambdaPhiseq=VARXlPhiseq, gran1Phi=VARXPhigran1, gran2Phi=VARXPhigran2,
                       lambdaBseq=VARXlBseq, gran1B=VARXBgran1, gran2B=VARXBgran2, eps=eps, max.iter=100,
                       T1.cutoff=cvcut, alpha=VARXalpha, type=VARXpen)

    VARXmodel <- HVARX_NEW_export_cpp(fullY=VARXdata$fullY, fullZ=VARXdata$fullZ, fullX=VARXdata$fullX,
                                      k=VARXdata$k, kX=VARXdata$kX, p=VARXdata$p, s=VARXdata$s,
                                      lambdaPhi=VARXcv$lPhi_oneSE, lambdaB=VARXcv$lB_oneSE,
                                      eps=eps, max_iter=100, alpha=VARXalpha, type=estim,
                                      Binit = matrix(0, VARXdata$k, VARXdata$kX*VARXdata$s),
                                      Phiinit =  matrix(0, VARXdata$k, VARXdata$k*VARXdata$p))
    Phihat <- VARXmodel$Phi
    phi0hat <- t(VARXmodel$phi0)
    Bhat <- VARXmodel$B
    out <- list("k"=k, "Y"=Y, "X"=X, "m"=m,"p"=p, "s"=s ,
                "Phihat"=Phihat, "Bhat"=Bhat, "phi0hat"=phi0hat,
                "exogenous_series_names"=exogenous_series_names,
                "endogenous_series_names"=endogenous_series_names,
                "lambdaPhi"=VARXcv$l1$lPhiseq, "lambdaB"=VARXcv$l1$lBseq,
                "lambdaPhi_opt"=VARXcv$lPhi_opt, "lambdaPhi_SEopt"=VARXcv$lPhi_oneSE,
                "lambdaB_opt"=VARXcv$lB_opt, "lambdaB_SEopt"=VARXcv$lB_oneSE,
                "MSFEcv"=VARXcv$MSFE_avg, "h"=h, selection = selection)

  }else { # Not time series cross validation
    # Grids regularization parameters
    if(is.null(VARXlPhiseq)){
      VARXlPhiseq <- PenaltyGrid(fullY=VARXdata$fullY, fullZ=VARXdata$fullZ, fullX=VARXdata$fullX, k=VARXdata$k, kX=VARXdata$kX, p=VARXdata$p, s=VARXdata$s,
                                  gran1Phi=VARXPhigran1, gran2Phi=VARXPhigran2, gran1B=VARXBgran1, gran2B=VARXBgran2, eps=eps,
                                  max.iter=100, alpha=VARXalpha, type=VARXpen)$lambdaPhiseq
    }
    if(is.null(VARXlBseq)){
      VARXlBseq <- PenaltyGrid(fullY=VARXdata$fullY, fullZ=VARXdata$fullZ, fullX=VARXdata$fullX, k=VARXdata$k, kX=VARXdata$kX, p=VARXdata$p, s=VARXdata$s,
                                gran1Phi=VARXPhigran1, gran2Phi=VARXPhigran2, gran1B=VARXBgran1, gran2B=VARXBgran2, eps=eps,
                                max.iter=100, alpha=VARXalpha, type=VARXpen)$lambdaBseq
    }
    l1 <- list(lPhiseq=sort(rep(VARXlPhiseq,length(VARXlBseq))), lBseq=rep(VARXlBseq,length(VARXlPhiseq)))
    l1.mat <- do.call(cbind, l1)

    Phihat <- array(NA, dim = c(VARXdata$k, VARXdata$k * VARXdata$p, length(l1$lPhiseq)))
    phi0hat <- array(NA, dim = c(VARXdata$k, 1, length(l1$lPhiseq)))
    Bhat <- array(NA, dim = c(VARXdata$k, VARXdata$kX * VARXdata$s, length(l1$lPhiseq)))
    for (i in 1:nrow(l1.mat)){
      VARXmodel <- HVARX_NEW_export_cpp(fullY=VARXdata$fullY, fullZ=VARXdata$fullZ, fullX=VARXdata$fullX,
                                       k=VARXdata$k, kX=VARXdata$kX, p=VARXdata$p, s=VARXdata$s,
                                       lambdaPhi=l1.mat[i, 1], lambdaB=l1.mat[i, 2],
                                       eps=eps, max_iter=100, alpha=VARXalpha, type=estim,
                                       Binit = matrix(0, VARXdata$k, VARXdata$kX*VARXdata$s),
                                       Phiinit =  matrix(0, VARXdata$k, VARXdata$k*VARXdata$p))
      Phihat[, , i] <- VARXmodel$Phi
      phi0hat[, , i] <- t(VARXmodel$phi0)
      Bhat[, , i] <- VARXmodel$B

      out <- list("k"=k, "Y"=Y, "X"=X, "m"=m,"p"=p, "s"=s ,
                  "Phihat"=Phihat, "Bhat"=Bhat, "phi0hat"=phi0hat,
                  "exogenous_series_names"=exogenous_series_names,
                  "endogenous_series_names"=endogenous_series_names,
                  "lambdaPhi"=l1$lPhiseq, "lambdaB"=l1$lBseq,
                  "lambdaPhi_opt"=NA, "lambdaPhi_SEopt"=NA,
                  "lambdaB_opt"=NA, "lambdaB_SEopt"=NA,
                  "MSFEcv"=NA, "h"=h, selection = selection)
    }
  }

  class(out) <- "bigtime.VARX"
  if (selection %in% c("bic", "aic", "hq")) out <- ic_selection(out, ic = selection, verbose = verbose)
  out
}

HVARXmodel <- function(Y, X, p, s, h=1, Yhat=NULL, cvY="default"){

  # Preliminaries
  k <- ncol(Y) # Number of Endogenous variables
  kX <- ncol(X) # Number of Exogenous variables
  cvYdata <- NULL

  # Lagged predictor matrices
  m <- max(s, p)
  DATAY <- embed(Y, dimension=m+h)
  if(!is.null(Yhat)){
    DATAYhat <- embed(Yhat, dimension=m+h)
    fullY <- DATAYhat[, 1:k] # in its columns Y1, Y2, ... Yq
  }else{
    fullY <- DATAY[, 1:k] # in its columns Y1, Y2, ... Yq
  }


  fullZNEW <- as.matrix(as.matrix(DATAY[, -c(1:k)])[, (1:((p+h-1)*k))])
  fullZ <- t(fullZNEW[, (ncol(fullZNEW)- k*p+1): ncol(fullZNEW)])

  DATAX <- embed(X, dimension=m+h)
  fullXNEW <- as.matrix(as.matrix(DATAX[, -c(1:kX)])[, (1: ((s+h-1)*kX))])
  fullX <- t(fullXNEW[, (ncol(fullXNEW)- kX*s+1): ncol(fullXNEW)])

  if(cvY=="predictY"){
    cvYdata <- DATAY[,1:k]
  }

  # Output
  out <- list("fullY"=fullY, "fullX"=fullX, "fullZ"=fullZ, "k"=k,"kX"=kX,"p"=p,"s"=s,
            "cvYdata"=cvYdata)
}

HVARX_cv<- function(Y, X, p, s, h=1, lambdaPhiseq=NULL, gran1Phi=10^2, gran2Phi=10, lambdaBseq=NULL, gran1B=10^2, gran2B=10, eps=10^-5, max.iter=100,
                    T1.cutoff=0.9, alpha=0, type="HLag"){

  # Get response and predictor matrices
  HVARXmodelFIT <- HVARXmodel(Y=Y, X=X, p=p, s=s, h=h, Yhat=NULL, cvY="default")
  k <- HVARXmodelFIT$k
  kX <- HVARXmodelFIT$kX
  fullY <- HVARXmodelFIT$fullY
  if (k==1){
    fullY <- matrix(fullY, ncol=1)
  }
  fullX <- HVARXmodelFIT$fullX
  fullZ <- HVARXmodelFIT$fullZ
  fullYtest <- HVARXmodelFIT$cvYdata

  # Grids regularization parameters
  if(is.null(lambdaPhiseq)){
    lambdaPhiseq <- PenaltyGrid(fullY=fullY, fullZ=fullZ, fullX=fullX, k=k, kX=kX, p=p, s=s,
                                    gran1Phi=gran1Phi, gran2Phi=gran2Phi, gran1B=gran1B, gran2B=gran2B, eps=eps,
                                    max.iter=max.iter, alpha=alpha, type=type)$lambdaPhiseq
  }
  if(is.null(lambdaBseq)){
    lambdaBseq <- PenaltyGrid(fullY=fullY, fullZ=fullZ, fullX=fullX, k=k, kX=kX, p=p, s=s,
                                  gran1Phi=gran1Phi, gran2Phi=gran2Phi, gran1B=gran1B, gran2B=gran2B, eps=eps,
                                  max.iter=max.iter, alpha=alpha, type=type)$lambdaBseq
  }
  l1 <- list(lPhiseq=sort(rep(lambdaPhiseq,length(lambdaBseq))), lBseq=rep(lambdaBseq,length(lambdaPhiseq)))

  # Time series cross-validation
  n <- nrow(fullY)
  T1 <- floor(T1.cutoff*n)
  tseq <- T1:(n-1)

  estim <- (type=="HLag")*2 + (type=="L1")*1
  my_cv <- HVARX_cvaux_cpp_loop(Y = fullY, Z = fullZ, X = fullX,  tseq = tseq, lambdaPhiseq = l1$lPhiseq , lambdaBseq = l1$lBseq, eps = eps,  max_iter = max.iter,
                                k = k,  kX = kX, p = p, s = s,  alpha = alpha,  estim = estim)

  MSFEmatrix <- my_cv$MSFE
  sparsitymatrix <- my_cv$sparsity

  MSFE_avg <- apply(MSFEmatrix, 2, mean)
  lPhiopt <- l1$lPhiseq[which.min(MSFE_avg)]
  lBopt <- l1$lBseq[which.min(MSFE_avg)]

  if(length(lPhiopt)==0){
    lPhiopt <- l1$lPhiseq[round(median(1:length(l1$lPhiseq)))]
  }else{
    if(is.na(lPhiopt)){
      lPhiopt <- l1$lPhiseq[round(median(1:length(l1$lPhiseq)))]
    }
  }
  if(length(lBopt)==0){
    lBopt <- l1$lBseq[round(median(1:length(l1$lBseq)))]
  }else{
    if(is.na(lBopt)){
      lBopt <- l1$lBseq[round(median(1:length(l1$lBseq)))]
    }
  }

  # One-standard error rule to determine optimal lambda values
  lPhi_oneSE <- NA
  lB_oneSE <- NA
  lphioptinit <- lPhiopt
  lBoptinit <- lBopt
  MSFE_sd <- apply(MSFEmatrix,2,sd)/sqrt(nrow(MSFEmatrix))
  MSFE_flagOK <- MSFE_avg<min(MSFE_avg,na.rm=T) + MSFE_sd[which.min(MSFE_avg)]
  MSFEnew <- MSFE_avg
  MSFEnew[!MSFE_flagOK] <- NA
  sparsitynew <- apply(sparsitymatrix,2,mean)
  sparsitynew[!MSFE_flagOK] <- NA
  sparsitynew[l1$lPhiseq<lPhiopt] <- NA
  sparsitynew[l1$lBseq<lBopt] <- NA

  lPhi_oneSE <- l1$lPhiseq[which.min(sparsitynew)]
  lB_oneSE <- l1$lBseq[which.min(sparsitynew)]

  if(length(lPhi_oneSE)==0){
    lPhi_oneSE <- lphioptinit
  }else{
    if(is.na(lPhi_oneSE)){
      lPhi_oneSE <- lphioptinit
    }
  }

  if(length(lB_oneSE)==0){
    lB_oneSE <- lBoptinit
  }else{
    if(is.na(lB_oneSE)){
      lB_oneSE <- lBoptinit
    }
  }

  flagPhi_oneSE <- lPhi_oneSE==min(lambdaPhiseq)
  flagB_oneSE <- lB_oneSE==min(lambdaBseq)
  if(lPhi_oneSE==min(lambdaPhiseq)){
    warning("Lower bound of lambdaPhi grid is selected: higher value of first granularity parameter is recommended as an input")
  }
  if(lB_oneSE==min(lambdaBseq)){
    warning("Lower bound of lambdaTheta/lambdaB grid is selected: higher value of first granularity parameter is recommended as an input")
  }

  out<-list("l1"=l1, "MSFEmatrix"=MSFEmatrix, "MSFE_avg"=MSFE_avg, "lPhi_oneSE"=lPhi_oneSE,
            "lB_oneSE"=lB_oneSE, "flag_oneSE"=c(flagPhi_oneSE,flagB_oneSE),
            "lPhi_opt"=lPhiopt, "lB_opt"=lBopt)
}

PenaltyGrid<- function(fullY, fullZ, fullX, k, kX, p, s, gran1Phi, gran2Phi, gran1B, gran2B, eps=eps,
                       max.iter=max.iter, alpha, type="HLag"){
  lambdamaxZ <- c()
  for (i in 1:k) {
    lambdamaxZ[i]  <- norm2(fullZ %*% fullY[, i])
  }
  lambdastartZ <- max(lambdamaxZ)

  lambdamaxX <- c()
  for (i in 1:k) {
    lambdamaxX[i]  <- norm2(fullX %*% fullY[, i])
  }
  lambdastartX <- max(lambdamaxX)


  SparsityLSearch <- SparsityLineSearch(fullY=fullY, fullZ=fullZ, fullX=fullX, k=k, kX=kX, p=p, s=s,
                                            lstartX=lambdastartX, lstartZ=lambdastartZ, eps=eps,
                                            max.iter=max.iter, alpha=alpha, type=type)

  lambdaPhiseq <- exp(seq(from = log(SparsityLSearch$lambdahZ), to = log(SparsityLSearch$lambdahZ/gran1Phi),
                          length = gran2Phi))
  lambdaBseq <- exp(seq(from = log(SparsityLSearch$lambdahX), to = log(SparsityLSearch$lambdahX/gran1B),
                        length = gran2B))
  out<-list("lambdaPhiseq"=lambdaPhiseq, "lambdaBseq"=lambdaBseq)
}

SparsityLineSearch<- function(fullY, fullZ, fullX, k, kX, p, s, lstartX, lstartZ, eps, max.iter, alpha, type="HLag"){
  lambdahX <- lstartX
  lambdalX <- 0
  lambdahZ <- lstartZ
  lambdalZ <- 0
  thresh <- 10

  estim <- (type=="HLag")*2 + (type=="L1")*1

  while(thresh>.00001)
  {
    lambdaX <- (lambdahX+lambdalX)/2
    lambdaZ <- (lambdahZ+lambdalZ)/2


    HVARXFIT <- HVARX_NEW_export_cpp(fullY=fullY, fullZ=fullZ, fullX=fullX, k=k, kX=kX, p=p, s=s,
                                     lambdaPhi=lambdaZ, lambdaB=lambdaX, eps=eps, max_iter=max.iter, alpha=alpha, type=estim, Binit = matrix(0, k, kX*s), Phiinit =  matrix(0, k, k*p))

    paramZ <- HVARXFIT$Phi
    paramX <- HVARXFIT$B

    if(max(abs(paramZ))==0)
    {
      lambdahZ <- lambdaZ

    }else{
      lambdalZ <- lambdaZ
    }

    if(max(abs(paramX))==0)
    {
      lambdahX <- lambdaX

    }else{
      lambdalX <- lambdaX
    }

    threshX <- max(abs(lambdahX-lambdalX))
    threshZ <- max(abs(lambdahZ-lambdalZ))
    thresh <- max(threshX,threshZ)

  }
  out <- list("lambdahX"=lambdahX, "lambdahZ"=lambdahZ)
}

