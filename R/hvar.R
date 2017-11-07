#' Sparse estimation of the Vector AutoRegressive (VAR) model
#' @param Y A \eqn{T} by \eqn{k} matrix of time series. If k=1, a univariate autoregressive model is estimated.
#' @param p User-specified maximum autoregressive lag order of the VAR. Typical usage is to have the program compute its own maximum lag order based on the time series length.
#' @param h Desired forecast horizon in time-series cross-validation procedure.
#' @param VARlseq User-specified grid of values for regularization parameter corresponding to sparse penalty. Typical usage is to have the program compute
#' its own grid. Supplying a grid of values overrides this. WARNING: use with care.
#' @param VARgran User-specified vector of granularity specifications for the penalty parameter grid:  First element specifies
#' how deep the grid should be constructed. Second element specifies how many values the grid should contain.
#' @param cvcut Proportion of observations used for model estimation in the time series cross-validation procedure. The remainder is used for forecast evaluation.
#' @param eps a small positive numeric value giving the tolerance for convergence in the proximal gradient algorithm.
#' @param VARalpha a small positive regularization parameter value corresponding to squared Frobenius penalty. The default is zero.
#' @param VARpen "HLag" (hierarchical sparse penalty) or "L1" (standard lasso penalty) penalization.
#' @export
#' @return A list with the following components
#' \item{Y}{\eqn{T} by \eqn{k} matrix of time series.}
#' \item{k}{Number of time series.}
#' \item{p}{Maximum autoregressive lag order of the VAR.}
#' \item{Phihat}{Matrix of estimated autoregressive coefficients of the VAR.}
#' \item{phi0hat}{vector of VAR intercepts.}
#' @references Nicholson William B., Bien Jacob and Matteson David S. (2017), "High Dimensional Forecasting via Interpretable Vector Autoregression"
#' arXiv preprint arXiv:1412.5250v2.
#' @examples
#' data(Y)
#' varfit <- sparsevar(Y) # sparse VAR
#' Y1 <- matrix(Y[,1], ncol=1)
#' arfit <- sparsevar(Y1) # sparse AR
sparsevar <- function(Y, p=NULL, VARpen="HLag", VARlseq=NULL, VARgran=NULL, VARalpha=0,
                      cvcut=0.9, h=1,  eps=1e-3){

  # Check Inputs
  if(!is.matrix(Y)){
    stop("Y needs to be a matrix of dimension T by k")
  }

  if(h<=0){
    stop("The forecast horizon h needs to be a strictly positive integer")
  }

  if(!is.null(p)){
    if(p<=0){
      stop("The maximum autoregressive order of VAR needs to be a strictly positive integer")
    }
  }

  if((!is.vector(VARlseq) & !is.null(VARlseq)) | length(VARlseq)==1){
    stop("The regularization parameter VARlseq needs to be a vector of length>1 or NULL otherwise")
  }

  if(any((VARgran<=0)==T)){
    stop("The granularity parameters needs to be a strictly positive integers")
  }

  if(!((cvcut<1) & (cvcut>0))){
    stop("cvcut needs to be a number between 0 and 1")
  }

  if(VARalpha<0){
    stop("The regularization paramter VARalpha needs to be a equal to zero or a small positive number")
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

  # Get optimal sparsity parameter via time series cross-validation
  VARcv <- HVAR_cv(Y=Y, p=p, h=h, lambdaPhiseq=VARlseq, gran1=VARgran1, gran2=VARgran2, T1.cutoff=cvcut,
                 eps=eps , VAR.alpha=VARalpha, type=VARpen)

  # Var estimation with selected regularization parameter
  VARdata <- HVARmodel(Y=Y, p=p, h=h)
  VARmodel <- HVAR(fullY=VARdata$fullY, fullZ=VARdata$fullZ, p=VARdata$p, k=VARdata$k,
             lambdaPhi=VARcv$lambda_opt_oneSE, eps=eps, VAR.alpha=VARalpha, type=VARpen)

  k <- ncol(Y)

  out <- list("k"=k, "Y"=Y, "p"=p, "Phihat"=VARmodel$Phi, "phi0hat"=VARmodel$phi)
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


HVAR<-function(fullY, fullZ, p, k, lambdaPhi, eps=1e-5, VAR.alpha=0, type="HLag"){

  # HVAR estimates
  if(type=="HLag"){
  if(k==1){
    fullY <- matrix(fullY, ncol=1)
  }
  Phi <- .HVARElemAlg(array(0, dim=c(k,k*p+1,1)), fullY, fullZ,
                      lambdaPhi, eps=eps, p, F, rep(1,p))
  }
  # L1 estimates
  if(type=="L1"){
    Phi <- .lassoVARFist(array(0,dim=c(k,k*p+1,1)), Z=fullZ, Y=fullY, gamm=lambdaPhi, eps=eps, p, F, rep(1,p))
  }
  Phi <- (1/(1 + VAR.alpha))*Phi
  resids <- t(t(fullY) - Phi[,,1]%*%rbind(rep(1, ncol(fullZ)), fullZ))
  Yhat <- t(Phi[,, 1]%*%rbind(rep(1, ncol(fullZ)), fullZ))

  # Output
  out <- list("Phi"=Phi[,-1,1], "phi0"=Phi[,1,1], "resids"=resids, "Yhat"=Yhat, "Y"=fullY)

}

HVAR_cv<-function(Y, p, h=1, lambdaPhiseq=NULL, gran1=20, gran2=10, T1.cutoff=0.9,
                  eps=1e-5 , VAR.alpha=0, type="HLag"){

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

  if(type=="HLag"){
  MSFEcv <-lapply(tseq, HVAR_cvaux, Y=fullY, Z=fullZ, gamm=lambdaPhiseq, eps=eps, p=p,
                  MN=F, C=rep(1,p), estim="HVARELEM", VAR.alpha=VAR.alpha)
  }

  if(type=="L1"){
    MSFEcv<-lapply(tseq, HVAR_cvaux, Y=fullY, Z=fullZ, gamm=lambdaPhiseq, eps=eps, p=p,
                   MN=F, C=rep(1,p), estim="Basic", VAR.alpha=VAR.alpha)
  }

  allresults <- c()
  for(i in 1:length(MSFEcv)){
    allresults <- rbind(allresults, MSFEcv[[i]]$MSFEs, MSFEcv[[i]]$sparsity)
  }

  MSFEmatrix <- allresults[seq(from=1, by=2, length=length(tseq)),]
  sparsitymatrix <- allresults[seq(from=2, by=2, length=length(tseq)),]
  if(length(tseq)==1){
    MSFEmatrix <- matrix(MSFEmatrix, nrow=1)
    sparsitymatrix <- matrix(sparsitymatrix, nrow=1)
    MSFE_avg <- apply(MSFEmatrix, 2, mean)
  }else{
    MSFE_avg <- apply(MSFEmatrix, 2, mean)
  }


  lambda_opt <- lambdaPhiseq[which.min(MSFE_avg)]

  if(length(lambda_opt)==0){
    lambda_opt <- median(lambdaPhiseq)
  }else{
    if(is.na(lambda_opt)){
      lambda_opt <- median(lambdaPhiseq)
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
  out <- list("lambda"=lambdaPhiseq, "lambda_opt_oneSE"=lambda_opt_oneSE,
            "MSFE_all"=MSFEmatrix, "MSFE_avg"=MSFE_avg, "flag_oneSE"=gridflag_oneSE)
}


HVAR_cvaux<-function(Y, Z, t, gamm, eps, p, MN, C, estim="HVARELEM", VAR.alpha){
  # Auxiliary function of HVAR_cv to select optimal sparsity parameter based on time-series cross-validation

  # Preliminaries
  k <- ncol(Y) # number of time series



  # Training sets
  trainY <- Y[(1:t), ]
  trainZ <- Z[, (1:t)]
  if(!is.matrix(trainZ)){trainZ <- matrix(trainZ, ncol=1)}
  if(!is.matrix(trainY)){trainY <- matrix(trainY, ncol=1)}

  # Prediction sets
  testY <- matrix(Y[t+1,], nrow=1)
  testZ <- matrix(Z[,t+1], ncol=1)


  # HVAR estimation
  Phi <- array(0,dim=c(k, k*p+1, length(gamm)))



  if(estim=="HVARELEM"){

  Phi <- .HVARElemAlg(Phi, trainY, trainZ, gamm, eps=eps, p, MN, C)


  }

  if(estim=="Basic"){
    Phi <- .lassoVARFist(Phi, Z=trainZ, Y=trainY, gamm=gamm, eps=eps,p,MN,C)


  }
  Phi <- (1/(1+VAR.alpha))*Phi



  MSFEs <- apply(Phi, 3, function(U,Y,Z){mean((t(t(Y) - U%*%rbind(rep(1, ncol(Z)), Z)))^2)}, Y=testY, Z=testZ)
  sparsity <- apply(Phi, 3, function(U){length(which(U[,-1]!=0))})

  if(!(k==1 & p==1)){
    CHECK_ALLZERO <- apply(Phi, 3, function(U){ length(which((abs(U[,-1]))!=0)) < dim(U)[1]})
    ALLZERO <- CHECK_ALLZERO
    sparsity[ALLZERO] <- NA
    MSFEs[ALLZERO] <- NA
  }


  out <- list("MSFEs"=MSFEs, "sparsity"=sparsity)
}

