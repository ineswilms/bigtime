#' Sparse Estimation of the Vector AutoRegressive Moving Average (VARMA) Model
#' @param Y A \eqn{T} by \eqn{k} matrix of time series. If k=1, a univariate autoregressive moving average model is estimated.
#' @param U A \eqn{T} by \eqn{k} matrix of (approximated) error terms. Typical usage is to have the program estimate a high-order VAR model (Phase I) to get approximated error terms U.
#' @param VARp User-specified maximum  autoregressive lag order of the PhaseI VAR. Typical usage is to have the program compute its own maximum lag order based on the time series length.
#' @param VARgran User-specified vector of granularity specifications for the penalty parameter grid of the PhaseI VAR:  First element specifies
#' how deep the grid should be constructed. Second element specifies how many values the grid should contain.
#' @param VARlseq User-specified grid of values for regularization parameter in the PhaseI VAR. Typical usage is to have the program compute
#' its own grid. Supplying a grid of values overrides this. WARNING: use with care.
#' @param VARpen "HLag" (hierarchical sparse penalty) or "L1" (standard lasso penalty) penalization in PhaseI VAR.
#' @param VARMAp User-specified maximum autoregressive lag order of the VARMA. Typical usage is to have the program compute its own maximum lag order based on the time series length.
#' @param VARMAq User-specified maximum moving average lag order of the VARMA. Typical usage is to have the program compute its own maximum lag order based on the time series length.
#' @param VARMAlPhiseq User-specified grid of values for regularization parameter corresponding to the autoregressive coefficients in the VARMA. Typical usage is to have the program compute
#' its own grid. Supplying a grid of values overrides this. WARNING: use with care.
#' @param VARMAPhigran User-specified vector of granularity specifications for the penalty parameter grid corresponding to the autoregressive coefficients in the VARMA:  First element specifies
#' how deep the grid should be constructed. Second element specifies how many values the grid should contain.
#' @param VARMAlThetaseq User-specified grid of values for regularization parameter corresponding to the moving average coefficients in the VARMA. Typical usage is to have the program compute
#' its own grid. Supplying a grid of values overrides this. WARNING: use with care.
#' @param VARMAThetagran User-specified vector of granularity specifications for the penalty parameter grid corresponding to the moving average coefficients in the VARMA:  First element specifies
#' how deep the grid should be constructed. Second element specifies how many values the grid should contain.
#' @param VARMApen "HLag" (hierarchical sparse penalty) or "L1" (standard lasso penalty) penalization in the VARMA.
#' @param VARMAalpha a small positive regularization parameter value corresponding to squared Frobenius penalty in  VARMA. The default is zero.
#' @param eps a small positive numeric value giving the tolerance for convergence in the proximal gradient algorithms.
#' @param h Desired forecast horizon in time-series cross-validation procedure.
#' @param cvcut Proportion of observations used for model estimation in the time series cross-validation procedure. The remainder is used for forecast evaluation.
#' @export
#' @return A list with the following components

#' \item{Y}{\eqn{T} by \eqn{k} matrix of time series.}
#' \item{U}{Matrix of (approximated) error terms.}
#' \item{k}{Number of time series.}
#' \item{VARp}{Maximum autoregressive lag order of the PhaseI VAR.}
#' \item{VARPhihat}{Matrix of estimated autoregressive coefficients of the Phase I VAR.}
#' \item{VARphi0hat}{Vector of Phase I VAR intercepts.}
#' \item{VARMAp}{Maximum autoregressive lag order of the VARMA.}
#' \item{VARMAq}{Maximum moving average lag order of the VARMA.}
#' \item{Phihat}{Matrix of estimated autoregressive coefficients of the VARMA.}
#' \item{Thetahat}{Matrix of estimated moving average coefficients of the VARMA.}
#' \item{phi0hat}{Vector of VARMA intercepts.}
#' \item{series_names}{names of time series}
#' \item{PhaseI_lambas}{Phase I sparsity parameter grid}
#' \item{PhaseI_MSFEcv}{MSFE cross-validation scores for each value of the sparsity parameter in the considered grid}
#' \item{PhaseI_lambda_opt}{Phase I Optimal value of the sparsity parameter as selected by the time-series cross-validation procedure}
#' \item{PhaseI_lambda_SEopt}{Phase I Optimal value of the sparsity parameter as selected by the time-series cross-validation procedure and after applying the one-standard-error rule}
#' \item{PhaseII_lambdaPhi}{Phase II sparsity parameter grid corresponding to Phi parameters}
#' \item{PhaseII_lambdaTheta}{Phase II sparsity parameter grid corresponding to Theta parameters}
#' \item{PhaseII_lambdaPhi_opt}{Phase II Optimal value of the sparsity parameter (corresponding to Phi parameters) as selected by the time-series cross-validation procedure}
#' \item{PhaseII_lambdaPhi_SEopt}{Phase II Optimal value of the sparsity parameter (corresponding to Theta parameters) as selected by the time-series cross-validation procedure and after applying the one-standard-error rule}
#' \item{PhaseII_lambdaTheta_opt}{Phase II Optimal value of the sparsity parameter (corresponding to Phi parameters) as selected by the time-series cross-validation procedure}
#' \item{PhaseII_lambdaTheta_SEopt}{Phase II Optimal value of the sparsity parameter (corresponding to Theta parameters) as selected by the time-series cross-validation procedure and after applying the one-standard-error rule}
#' \item{PhaseII_MSFEcv}{Phase II MSFE cross-validation scores for each value in the two-dimensional sparsity grid}
#' \item{h}{Forecast horizon h}
#' @references Wilms Ines, Sumanta Basu, Bien Jacob and Matteson David S. (2017), "Sparse Identification and Estimation of High-Dimensional Vector AutoRegressive Moving Averages"
#' arXiv preprint <arXiv:1707.09208>.
#' @seealso \link{lagmatrix} and \link{directforecast}
#' @examples
#' data(Y)
#' VARMAfit <- sparseVARMA(Y) # sparse VARMA
#' y <- matrix(Y[,1], ncol=1)
#' ARMAfit <- sparseVARMA(y) # sparse ARMA
sparseVARMA <- function(Y, U=NULL,  VARp=NULL, VARpen="HLag", VARlseq=NULL, VARgran=NULL,
                        VARMAp=NULL, VARMAq=NULL, VARMApen="HLag",VARMAlPhiseq=NULL, VARMAPhigran=NULL,
                        VARMAlThetaseq=NULL, VARMAThetagran=NULL, VARMAalpha=0,
                        h=1, cvcut=0.9, eps=10^-3){

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

  if(!is.matrix(U) & !is.null(U)){

    if(is.vector(U) & length(U)>1){
      U <- matrix(U, ncol=1)
    }else{
      stop("U needs to be a matrix of dimension T by k or NULL otherwise.")
    }

  }

  if(!is.null(U)){
    if(nrow(U)!=nrow(Y)){
      stop("U and Y need to have the same number of rows (observations). Otherwise leave U NULL.")
    }
  }

  if(!is.null(colnames(Y))){ # time series names
    series_names <- colnames(Y)
  } else {
    series_names <- NULL
  }

  if(h<=0){
    stop("The forecast horizon h needs to be a strictly positive integer")
  }

  if(!is.null(VARp)){
    if(VARp<=0){
      stop("The maximum autoregressive order of the PhaseI VAR needs to be a strictly positive integer")
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

  if(!is.null(VARMAp)){
    if(VARMAp<=0){
      stop("The maximum autoregressive order of the VARMA needs to be a strictly positive integer")
    }
  }

  if(!is.null(VARMAq)){
    if(VARMAq<=0){
      stop("The maximum moving average order of the VARMA needs to be a strictly positive integer")
    }
  }

  if( (!is.vector(VARMAlPhiseq) & !is.null(VARMAlPhiseq)) | length(VARMAlPhiseq)==1){
    stop("The regularization parameter grid VARMAlPhiseq needs to be a vector of length > 1 or NULL otherwise")
  }

  if(any((VARMAPhigran<=0)==T)){
    stop("The granularity parameters need to be a strictly positive integer")
  }


  if((!is.vector(VARMAlThetaseq) & !is.null(VARMAlThetaseq)) | length(VARMAlThetaseq)==1){
    stop("The regularization parameter VARMAlThetaseq needs to be a vector of length >1 or NULL otherwise")
  }

  if(any((VARMAThetagran<=0)==T)){
    stop("The granularity parameters need to be a strictly positive integer")
  }

  if(VARMAalpha<0){
    stop("The regularization paramter VARMAalpha needs to be a equal to zero or a small positive number")
  }

  # if(VARalpha<0){
  #   stop("The regularization paramter VARalpha needs to be a equal to zero or a small positive number")
  # }

  if(!is.element(VARpen, c("HLag", "L1"))){
    stop("The type of penalization VARpen needs to be either HLag (hierarchical sparse penalization) or L1 (standard lasso penalization)")
  }

  if(!is.null(VARp)){
    if(VARpen=="HLag" & VARp<=1){
      stop("HLag penalization in PhaseI VAR is only supported for p>1. Use L1 as VARpen instead")
    }
  }


  if(!is.element(VARMApen, c("HLag", "L1"))){
    stop("The type of penalization VARMApen needs to be either HLag (hierarchical sparse penalization) or L1 (standard lasso penalization)")
  }

  if(!is.null(VARMAq) | !is.null(VARMAp)){

    if(is.null(VARMAp)){
      VARMAp <- floor(1.5*sqrt(nrow(Y))/2)
    }

    if(is.null(VARMAq)){
      VARMAq <- floor(1.5*sqrt(nrow(Y))/2)
    }

    if(VARMApen=="HLag" & (VARMAp<=1 |VARMAq<=1)){
      stop("HLag penalization in VARMA is only supported for p and q larger than 1. Use L1 as VARMApen instead")
    }
  }


  if(eps<=0){
    stop("The convergence tolerance parameter eps needs to be a small positive number")
  }

  .check_if_standardised(Y)

  # Initialization
  HVARcvFIT  <- HVARFIT <-  NA
  HVARXcvFIT <- HVARXFIT <- NA
  VARPhi <- VARphi0 <- NA
  Phi <- Theta <- phi0 <- resid <- NA

  #### START PHASE I ####
  # Determine maximum order of VAR as function of sample size
  if(is.null(VARp)){
    VARp <- floor(1.5*sqrt(nrow(Y)))
  }
  VARp <- round(VARp)

  # Regularization grid
  if(is.null(VARgran)){
    VARgran1 <- 10^5
    VARgran2 <- 10
  }else{
    VARgran1 <- VARgran[1]
    VARgran2 <- VARgran[2]
  }

  if(!is.null(VARlseq)){
    VARgran1 <- max(VARlseq)/min(VARlseq)
    VARgran2 <- length(VARlseq)
  }

  if(!is.null(U)){
    HVARcvFIT <- list()
    HVARcvFIT$lambda <- NA
    HVARcvFIT$MSFE_avg<- NA
    HVARcvFIT$lambda_opt_oneSE <- NA
    HVARcvFIT$lambda_opt <- NA
  }

  # Fit a VAR(VARp) using HLag penalty to approximate the residuals
  if(is.null(U)){
    ### HVAR ESTIMATION

    HVARcvFIT <- HVAR_cv(Y = Y, p = VARp, h = h, lambdaPhiseq=VARlseq, gran1 = VARgran1, gran2=VARgran2, T1.cutoff=cvcut, eps=eps, type=VARpen)


    lambdaFLAG <- HVARcvFIT$flag
    HVARmodelFIT <- HVARmodel(Y=Y, p=VARp, h=h)
    if(HVARmodelFIT$k==1){
      HVARmodelFIT$fullY <- as.matrix(HVARmodelFIT$fullY)
    }

    HVARFIT <- HVAR(fullY=HVARmodelFIT$fullY, fullZ=HVARmodelFIT$fullZ, p=HVARmodelFIT$p, k=HVARmodelFIT$k, lambdaPhi=HVARcvFIT$lambda_opt_oneSE, eps=eps, type=VARpen)


    VARPhi <- HVARFIT$Phi
    VARphi0 <- HVARFIT$phi0

    # Obtain residuals
    HVARresid <- HVARFIT$resids
    nrowY <- nrow(Y)
    nrowresid <- nrow(HVARresid)
    nrowdiff <- nrowY - nrowresid
    U <- matrix(NA, ncol=ncol(Y), nrow=nrow(Y))
    U[-c(1:nrowdiff),] <- HVARresid
    U[1:nrowdiff,] <- HVARresid[1:nrowdiff,]
    U <- scale(U, center=T, scale=F)

  }else{
    HVARFIT <- NULL
  }

  #### START PHASE II ####
  ### Determine maximum AR and MA order of VARMA as a function of time series length T
  if(is.null(VARMAp)){
    VARMAp <- floor(1.5*sqrt(nrow(Y))/2)
  }

  if(is.null(VARMAq)){
    VARMAq <- floor(1.5*sqrt(nrow(Y))/2)
  }

  VARMAp <- round(VARMAp)
  VARMAq <- round(VARMAq)

  # Regularization grid
  if(is.null(VARMAPhigran)){
    VARMAPhigran1 <- 10^2
    VARMAPhigran2 <- 10
  }else{
    VARMAPhigran1 <- VARMAPhigran[1]
    VARMAPhigran2 <- VARMAPhigran[2]
  }

  if(!is.null(VARMAlPhiseq)){
    VARMAPhigran1 <- max(VARMAlPhiseq)/min(VARMAlPhiseq)
    VARMAPhigran2 <- length(VARMAlPhiseq)
  }

  if(is.null(VARMAThetagran)){
    VARMAThetagran1 <- 10^2
    VARMAThetagran2 <- 10
  }else{
    VARMAThetagran1 <- VARMAThetagran[1]
    VARMAThetagran2 <- VARMAThetagran[2]
  }

  if(!is.null(VARMAlThetaseq)){
    VARMAThetagran1 <- max(VARMAlThetaseq)/min(VARMAlThetaseq)
    VARMAThetagran2 <- length(VARMAlThetaseq)
  }



  ### Fit VARMA
  HVARXcvFIT <- HVARX_cv(Y=Y, X=U, p=VARMAp, s=VARMAq, h=h, lambdaPhiseq=VARMAlPhiseq, gran1Phi=VARMAPhigran1, gran2Phi=VARMAPhigran2,
                             lambdaBseq=VARMAlThetaseq, gran1B=VARMAThetagran1, gran2B=VARMAThetagran2, eps=eps, max.iter=100,
                             T1.cutoff=cvcut, alpha=VARMAalpha, type=VARMApen)

  HVARXmodelFIT <- HVARXmodel(Y=Y, X=U, p=VARMAp, s=VARMAq,h=h)

  estim <- (VARMApen=="HLag")*2 + (VARMApen=="L1")*1
  if(HVARXmodelFIT$k==1){
    HVARXmodelFIT$fullY <- as.matrix(HVARXmodelFIT$fullY)
  }
  HVARXFIT <- HVARX_NEW_export_cpp(fullY=HVARXmodelFIT$fullY, fullZ=HVARXmodelFIT$fullZ, fullX=HVARXmodelFIT$fullX,
                                   k=HVARXmodelFIT$k, kX=HVARXmodelFIT$kX, p=HVARXmodelFIT$p, s=HVARXmodelFIT$s,
                                   lambdaPhi=HVARXcvFIT$lPhi_oneSE, lambdaB=HVARXcvFIT$lB_oneSE,
                                   eps=eps, max_iter=100, alpha=VARMAalpha, type=estim, Binit = matrix(0, HVARXmodelFIT$k, HVARXmodelFIT$kX*HVARXmodelFIT$s), Phiinit =  matrix(0, HVARXmodelFIT$k, HVARXmodelFIT$k*HVARXmodelFIT$p))

  Phi <- HVARXFIT$Phi
  phi0 <- HVARXFIT$phi0
  Theta <- HVARXFIT$B
  resid <- HVARXFIT$resids
  U <- matrix(U, ncol=ncol(U), nrow=nrow(U))
  k <- ncol(Y)

  # Output
  out<-list("k"=k, "Y"=Y, "VARp"=VARp, "VARPhihat"=VARPhi, "VARphi0hat"=VARphi0, "U"=U,
            "VARMAp"=VARMAp, "VARMAq"=VARMAq,
            "Phihat"=Phi, "Thetahat"=Theta, "phi0hat"=t(phi0),
            "series_names"=series_names,
            "PhaseI_lambdas"=HVARcvFIT$lambda, "PhaseI_MSFEcv"=HVARcvFIT$MSFE_avg,
            "PhaseI_lambda_SEopt"=HVARcvFIT$lambda_opt_oneSE,"PhaseI_lambda_opt"=HVARcvFIT$lambda_opt,
            "PhaseII_lambdaPhi"=HVARXcvFIT$l1$lPhiseq, "PhaseII_lambdaTheta"=HVARXcvFIT$l1$lBseq,
            "PhaseII_lambdaPhi_opt"=HVARXcvFIT$lPhi_opt, "PhaseII_lambdaPhi_SEopt"=HVARXcvFIT$lPhi_oneSE,
            "PhaseII_lambdaTheta_opt"=HVARXcvFIT$lB_opt, "PhaseII_lambdaTheta_SEopt"=HVARXcvFIT$lB_oneSE,
            "PhaseII_MSFEcv"=HVARXcvFIT$MSFE_avg, "h"=h)
  class(out) <- "bigtime.VARMA"
  out
}
