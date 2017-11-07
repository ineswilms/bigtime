#' Sparse estimation of the Vector AutoRegressive  with exogenous variables X (VARX) model
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
#' @references Wilms Ines, Sumanta Basu, Bien Jacob and Matteson David S. (2017), "Interpretable Vector AutoRegressions with
#' Exogenous Time Series" arXiv preprint.
#' @seealso \link{lagmatrix} and \link{directforecast}
#' @examples
#' data(Y)
#' data(X)
#' varxfit <- sparsevarx(Y=Y, X=X) # sparse VARX
#' Y1 <- matrix(Y[,1], ncol=1)
#' arxfit <- sparsevarx(Y=Y1, X=X) # sparse ARX
sparsevarx <- function(Y, X, p=NULL, s=NULL, VARXpen="HLag", VARXlPhiseq=NULL, VARXPhigran=NULL,
                        VARXlBseq=NULL,  VARXBgran=NULL, VARXalpha=0, h=1, cvcut=0.9, eps=10^-3){

  # Check inputs
  # Check Inputs
  if(!is.matrix(Y)){

    if(is.vector(Y) & length(Y)>1){
      Y <- matrix(Y, ncol=1)
    }else{
      stop("Y needs to be a matrix of dimension T by k")
    }

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

  if( (!is.vector(VARXlPhiseq) & !is.null(VARXlPhiseq)) | length(VARXlPhiseq)==1){
    stop("The regularization parameter grid VARXlPhiseq needs to be a vector of length > 1 or NULL otherwise")
  }

  if(any((VARXPhigran<=0)==T)){
    stop("The granularity parameters need to be a strictly positive integer")
  }


  if((!is.vector(VARXlBseq) & !is.null(VARXlBseq)) | length(VARXlBseq)==1){
    stop("The regularization parameter VARXlBseq needs to be a vector of length >1 or NULL otherwise")
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

  if(!is.null(VARXPhigran)){
    VARXPhigran1 <- max(VARXPhigran)/min(VARXPhigran)
    VARXPhigran2 <- length(VARXPhigran)
  }

  if(is.null(VARXBgran)){
    VARXBgran1 <- 10^2
    VARXBgran2 <- 10
  }else{
    VARXBgran1 <- VARXBgran[1]
    VARXBgran2 <- VARXBgran[2]
  }

  if(!is.null(VARXBgran)){
    VARXBgran1 <- max(VARXBgran)/min(VARXBgran)
    VARXBgran2 <- length(VARXBgran)
  }

  # Get optimal sparsity parameter via time series cross-validation
  VARXcv <- HVARX_cv(Y=Y, X=X, p=p, s=s, h=h, lambdaPhiseq=VARXlPhiseq, gran1Phi=VARXPhigran1, gran2Phi=VARXPhigran2,
                     lambdaBseq=VARXlBseq, gran1B=VARXBgran1, gran2B=VARXBgran2, eps=eps, max.iter=100,
                     T1.cutoff=cvcut, alpha=VARXalpha, type=VARXpen)

  # VarX estimation with selected regularization parameter
  VARXdata <- HVARXmodel(Y=Y, X=X, p=p, s=s, h=h)
  VARXmodel <- HVARX(fullY=VARXdata$fullY, fullZ=VARXdata$fullZ, fullX=VARXdata$fullX,
                     k=VARXdata$k, kX=VARXdata$kX, p=VARXdata$p, s=VARXdata$s,
                     lambdaPhi=VARXcv$lPhi_oneSE, lambdaB=VARXcv$lB_oneSE,
                     eps=eps, max.iter=100, alpha=VARXalpha, type=VARXpen)

  k <- ncol(Y)
  m <- ncol(X)
  out <- list("k"=k, "Y"=Y, "X"=X, "m"=m,"p"=p, "s"=s ,
              "Phihat"=VARXmodel$Phi, "Bhat"=VARXmodel$B, "phi0hat"=VARXmodel$phi0)

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


HVARX <- function(fullY, fullZ, fullX, k, kX, p, s, lambdaPhi, lambdaB,
                  eps=10^-5, max.iter=100, alpha=0, type="HLag"){

  if(k==1){
    fullY <- matrix(fullY, ncol=1)
  }
  # Transformation to remove intercept
  YMean <- c(apply(fullY, 2, mean))
  ZMean <- c(apply(fullZ, 1, mean))
  XMean <- c(apply(fullX, 1, mean))
  Y <- fullY - matrix(c(rep(1, nrow(fullY))), ncol = 1) %*% matrix(c(apply(fullY, 2, mean)), nrow = 1)
  Z <- fullZ - c(apply(fullZ, 1, mean)) %*% t(c(rep(1, ncol(fullZ))))
  X <- fullX - c(apply(fullX, 1, mean)) %*% t(c(rep(1, ncol(fullX))))

  # Proximal gradient algorithm
  ZX <- rbind(Z,X)
  tk <- 1/max(Mod(eigen(ZX %*% t(ZX))$values))
  nbr.series <- c(1:k)
  if(type=="HLag"){
  PHIBETAprox <- lapply(nbr.series, FUN=proxHVARX, Ydata=Y, Zdata=Z, Xdata=X, tk=tk, lambdaPhi=lambdaPhi,
                      lambdaB=lambdaB, k=k, p=p, kX=kX, s=s, eps=eps, max.iter=max.iter, alpha=alpha)
  }

  if(type=="L1"){
    PHIBETAprox<-lapply(nbr.series, FUN=proxBasic, Ydata=Y, Zdata=Z, Xdata=X, tk=tk, lambdaPhi=lambdaPhi,
                        lambdaB=lambdaB, k=k, p=p, kX=kX, s=s, eps=eps, max.iter=max.iter, alpha=alpha)
  }

  PHIBETA <- matrix(unlist(PHIBETAprox), ncol=k*p+kX*s, nrow=k, byrow=T)
  PHI <- PHIBETA[, 1:(k*p)]
  B <- PHIBETA[, -c(1:(k*p))]

  phi0 <- YMean - PHI%*%ZMean - B%*%XMean
  if(k==1){
    resids <- t(t(fullY) - c(phi0,PHI)%*%rbind(rep(1, ncol(fullZ)), fullZ)- B%*%fullX)
  }else{
    resids <- t(t(fullY) - cbind(phi0,PHI)%*%rbind(rep(1, ncol(fullZ)), fullZ)- B%*%fullX)
  }


  # Output
  out <- list("Phi"=PHI, "B"=B, "phi0"=phi0, "resids"=resids,
            "fullY"=fullY, "fullZ"=fullZ, "fullX"=fullX)
}

proxHVARX <- function(i.series, Ydata, Zdata, Xdata, tk, lambdaPhi, lambdaB,
                      k, p, kX, s, eps=10^-5, max.iter=100, alpha=0){
  # Auxiliary Function: proximal gradient algorithm for HVARX estimation
  # Initialization
  it <- 3
  PhiOLD <- PhiOLDOLD <- matrix(0, ncol=k*p, nrow=1)
  BOLD <- BOLDOLD <- matrix(0, ncol=kX*s, nrow=1)
  thresh <- eps*10

  # Iterate until convergence
  while( (thresh>eps) & (it<max.iter)){
    phi <- PhiOLD + ((it-2)/(it+1))*(PhiOLD-PhiOLDOLD)
    beta <- BOLD + ((it-2)/(it+1))*(BOLD-BOLDOLD)
    gradientZ <- -(matrix(Ydata[,i.series], nrow=1)- matrix(phi,nrow=1)%*%Zdata - matrix(beta,nrow=1)%*%Xdata)%*%t(Zdata)
    gradientX <- -(matrix(Ydata[,i.series], nrow=1)- matrix(phi,nrow=1)%*%Zdata - matrix(beta,nrow=1)%*%Xdata)%*%t(Xdata)
    argZ <- phi - tk*gradientZ
    argX <- beta - tk*gradientX

    proxPhi <- prox2HVAR(v=argZ, lambda=tk*lambdaPhi, k=k, p=p)
    proxB <- prox2HVAR(v=argX, lambda=tk*lambdaB, k=kX, p=s)

    threshPhi <- max(abs(proxPhi - phi))
    threshB <- max(abs(proxB - beta))
    thresh <- max(threshPhi, threshB)

    PhiOLDOLD <- PhiOLD
    BOLDOLD <- BOLD
    PhiOLD <- proxPhi
    BOLD <- proxB
    it <- it+1
  }
  proxPhi <- (1/(1+alpha))*proxPhi
  proxB <- (1/(1+alpha))*proxB

  out <- list("Phi"=proxPhi, "B"=proxB)
}

proxBasic<-function(i.series,Ydata,Zdata,Xdata,tk,lambdaPhi,lambdaB,k,p,kX,s,eps=10^-5,max.iter=100, alpha=0){
  # Auxiliary Function: proximal gradient algorithm for l1-VARX estimation

  # Initialization
  it <- 3
  PhiOLD <- PhiOLDOLD <- matrix(0, ncol=k*p, nrow=1)
  BOLD <- BOLDOLD <- matrix(0, ncol=kX*s, nrow=1)
  thresh <- eps*10

  if(k==1 & p==1){
    Zdata <- matrix(Zdata, nrow=1)
  }

  if(kX==1 & s==1){
    Xdata <- matrix(Xdata, nrow=1)
  }
  # Iterate until convergence
  while( (thresh>eps) & (it<max.iter)){
    phi <- PhiOLD + ((it-2)/(it+1))*(PhiOLD-PhiOLDOLD)
    beta <- BOLD + ((it-2)/(it+1))*(BOLD-BOLDOLD)


    gradientZ <- -(matrix(Ydata[,i.series], nrow=1) - matrix(phi,nrow=1)%*%Zdata - matrix(beta,nrow=1)%*%Xdata)%*%t(Zdata)
    gradientX <- -(matrix(Ydata[,i.series], nrow=1) - matrix(phi,nrow=1)%*%Zdata - matrix(beta,nrow=1)%*%Xdata)%*%t(Xdata)
    argZ <- phi - tk*gradientZ
    argX <- beta - tk*gradientX

    proxPhi <- unlist(lapply(argZ, ST1a, gam=tk*lambdaPhi))
    proxB <- unlist(lapply(argX, ST1a, gam=tk*lambdaB))

    threshPhi <- max(abs(proxPhi - phi))
    threshB <- max(abs(proxB - beta))
    thresh <- max(threshPhi, threshB)

    PhiOLDOLD <- PhiOLD
    BOLDOLD <- BOLD
    PhiOLD <- proxPhi
    BOLD <- proxB
    it <- it+1
  }
  proxPhi <- (1/(1+alpha))*proxPhi
  proxB <- (1/(1+alpha))*proxB

  out<-list("Phi"=proxPhi, "B"=proxB)
}


HVARX_cv<-function(Y, X, p, s, h=1, lambdaPhiseq=NULL, gran1Phi=20, gran2Phi=10,
                   lambdaBseq=NULL, gran1B=20, gran2B=10, eps=10^-5, max.iter=100,
                   T1.cutoff=0.9, YhatPhaseI=NULL, cvYtype="default", alpha=0, type="HLag"){

  # Get response and predictor matrices
  HVARXmodelFIT <- HVARXmodel(Y=Y, X=X, p=p, s=s, h=h, Yhat=YhatPhaseI, cvY=cvYtype)
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
  MSFE_FIT <- lapply(tseq, HVARX_cvaux, fullY=fullY, fullX=fullX, fullZ=fullZ, eps=eps, max.iter=max.iter,
                     k=k, kX=kX, p=p, s=s, l1=l1, lambdaPhiseq=lambdaPhiseq, lambdaBseq=lambdaBseq,
                      cvY=cvYtype, fullYtest=fullYtest, alpha=alpha, type=type)
  allresults <- c()
  for(i in 1:length(MSFE_FIT)){
    allresults <- rbind(allresults, MSFE_FIT[[i]])
  }
  MSFEmatrix <- allresults[seq(from=1, by=3, length=length(tseq)),]
  sparsitymatrix <- allresults[seq(from=2, by=3, length=length(tseq)),]
  MSFEmatrix_relax <- allresults[seq(from=3, by=3, length=length(tseq)),]

  if(length(tseq)==1){
    MSFEmatrix <- matrix(MSFEmatrix, nrow=1)
    sparsitymatrix <- matrix(sparsitymatrix, nrow=1)
    MSFEmatrix_relax <- matrix(MSFEmatrix_relax, nrow=1)
    MSFE_avg <- apply(MSFEmatrix, 2, mean)
  }else{
    MSFE_avg <- apply(MSFEmatrix, 2, mean)
  }

  lPhiopt <- l1$lPhiseq[which.min(MSFE_avg)]
  lBopt <- l1$lBseq[which.min(MSFE_avg)]

  if(length(lPhiopt)==0){
    lPhiopt <- median(l1$lPhiseq)
  }else{
    if(is.na(lPhiopt)){
      lPhiopt <- median(l1$lPhiseq)
    }
  }
  if(length(lBopt)==0){
    lBopt <- median(l1$lBseq)
  }else{
    if(is.na(lBopt)){
      lBopt <- median(l1$lBseq)
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
            "lB_oneSE"=lB_oneSE, "flag_oneSE"=c(flagPhi_oneSE,flagB_oneSE))
}

HVARX_cvaux<-function(t, fullY, fullX, fullZ, eps, max.iter, k, kX, p, s, l1, lambdaPhiseq, lambdaBseq,
                      cvY="default", fullYtest=NULL, alpha, type="HLag"){
  # Auxiliary function of HVARX_cv to select optimal sparsity parameters based on time-series cross-validation

  # Training sets

  if(k==1){
    trainY <- fullY[(1:t),]
    trainY <- matrix(trainY, ncol=1)
  }else{
    trainY <- fullY[(1:t),]
  }

  if(s==1 & kX==1){
    trainX <- fullX[,(1:t)]
    trainX <- matrix(trainX, nrow=1)
  }else{
    trainX <- fullX[,(1:t)]
  }

  if(p==1 & k==1){
    trainZ <- fullZ[,(1:t)]
    trainZ <- matrix(trainZ, nrow=1)
  }else{
    trainZ <- fullZ[,(1:t)]
  }


  # Test sets
  testY <- matrix(fullY[t+1,], nrow=1)
  if(cvY=="predictY"){
    testY <- matrix(fullYtest[t,], nrow=1)
  }
  testX <- matrix(fullX[,t+1], ncol=1)
  testZ <- matrix(fullZ[,t+1], ncol=1)


  morearg <- list(Ytrain=trainY, Xtrain=trainX, Ztrain=trainZ, Ytest=testY, Xtest=testX, Ztest=testZ,
                  eps=eps, max.iter=max.iter, k=k, kX=kX, p=p, s=s, alpha=alpha, type=type)
  MSFEscv <- mapply(HVARX_MSFE, l1$lPhiseq, l1$lBseq, MoreArgs = morearg)
  allresults <- matrix(unlist(MSFEscv), ncol=length(l1$lPhiseq), nrow=3)
  return(allresults)
}

HVARX_MSFE<-function(Ytrain, Xtrain, Ztrain, Ytest, Xtest, Ztest, lamPhi, lamB,
                     eps, max.iter, k, kX, p, s, alpha, type="HLag"){
  # Auxiliary function to select optimal sparsity parameters based on time-series cross-validation

  HVARXFIT <- HVARX(fullY=Ytrain, fullZ=Ztrain, fullX=Xtrain, k=k, kX=kX, p=p, s=s,
                    lambdaPhi=lamPhi, lambdaB=lamB, eps=eps, max.iter=max.iter, alpha=alpha, type=type)
  Phi_FIT <- HVARXFIT$Phi
  Phi0_FIT <- HVARXFIT$phi0
  B_FIT <- HVARXFIT$B

  MSFErelax <- NA
  sparsity <- length(which(Phi_FIT!=0)) + length(which(B_FIT!=0))
  if(all(B_FIT==0)|all(Phi_FIT==0)){ # Check all zero solution -> exclude from options
    MSFE <- NA
    sparsity <- NA
  }else{

    if(k==1){
      MSFE <- mean(t( t(Ytest)- c(Phi0_FIT, Phi_FIT)%*%rbind(rep(1, ncol(Ztest)), Ztest) - B_FIT%*%Xtest)^2)
    }else{
      MSFE <- mean(t( t(Ytest)- cbind(Phi0_FIT, Phi_FIT)%*%rbind(rep(1, ncol(Ztest)), Ztest) - B_FIT%*%Xtest)^2)
    }

    MSFErelax <- NA
  }

  out<-list("MSFE"=MSFE, "sparsity"=sparsity, "MSFErelax"=MSFErelax)
}

PenaltyGrid<-function(fullY, fullZ, fullX, k, kX, p, s, gran1Phi, gran2Phi, gran1B, gran2B, eps=eps,
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


  SparsityLSearch<-SparsityLineSearch(fullY=fullY, fullZ=fullZ, fullX=fullX, k=k, kX=kX, p=p, s=s,
                                      lstartX=lambdastartX, lstartZ=lambdastartZ, eps=eps,
                                      max.iter=max.iter, alpha=alpha, type=type)

  lambdaPhiseq <- exp(seq(from = log(SparsityLSearch$lambdahZ), to = log(SparsityLSearch$lambdahZ/gran1Phi),
                          length = gran2Phi))
  lambdaBseq <- exp(seq(from = log(SparsityLSearch$lambdahX), to = log(SparsityLSearch$lambdahX/gran1B),
                        length = gran2B))
  out<-list("lambdaPhiseq"=lambdaPhiseq, "lambdaBseq"=lambdaBseq)
}

SparsityLineSearch<-function(fullY, fullZ, fullX, k, kX, p, s, lstartX, lstartZ, eps, max.iter, alpha, type="HLag"){
  lambdahX <- lstartX
  lambdalX <- 0
  lambdahZ <- lstartZ
  lambdalZ <- 0
  thresh <- 10
  while(thresh>.00001)
  {
    lambdaX <- (lambdahX+lambdalX)/2
    lambdaZ <- (lambdahZ+lambdalZ)/2

    HVARXFIT <- HVARX(fullY=fullY, fullZ=fullZ, fullX=fullX, k=k, kX=kX, p=p, s=s,
                      lambdaPhi=lambdaZ, lambdaB=lambdaX, eps=eps, max.iter=max.iter, alpha=alpha, type=type)
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

