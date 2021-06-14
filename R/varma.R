#' Sparse Estimation of the Vector AutoRegressive Moving Average (VARMA) Model
#' @param Y A \eqn{T} by \eqn{k} matrix of time series. If k=1, a univariate autoregressive moving average model is estimated.
#' @param U A \eqn{T} by \eqn{k} matrix of (approximated) error terms. Typical usage is to have the program estimate a high-order VAR model (Phase I) to get approximated error terms U.
#' @param VARp User-specified maximum  autoregressive lag order of the PhaseI VAR. Typical usage is to have the program compute its own maximum lag order based on the time series length.
#' @param VARgran User-specified vector of granularity specifications for the penalty parameter grid of the PhaseI VAR:  First element specifies
#' how deep the grid should be constructed. Second element specifies how many values the grid should contain.
#' @param VARlseq User-specified grid of values for regularization parameter in the PhaseI VAR. Typical usage is to have the program compute
#' its own grid. Supplying a grid of values overrides this. WARNING: use with care.
#' @param VARpen "HLag" (hierarchical sparse penalty) or "L1" (standard lasso penalty) penalization in PhaseI VAR.
#' @param VARselection Selection procedure for the first stage. Default is time series Cross-Validation. Alternatives are BIC, AIC, HQ
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
#' @param VARMAselection selection procedure in the second stage. Default is "none"; Alternatives are cv, bic, aic, hq
#' @param eps a small positive numeric value giving the tolerance for convergence in the proximal gradient algorithms.
#' @param h Desired forecast horizon in time-series cross-validation procedure.
#' @param cvcut Proportion of observations used for model estimation in the time series cross-validation procedure. The remainder is used for forecast evaluation.
#' @param check_std Check whether data is standardised. Default is TRUE and is not recommended to be changed
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
sparseVARMA <- function(Y, U=NULL,  VARp=NULL, VARpen="HLag", VARlseq=NULL,
                        VARgran=NULL, VARselection = c("cv", "bic", "aic", "hq"),
                        VARMAp=NULL, VARMAq=NULL, VARMApen="HLag",
                        VARMAlPhiseq=NULL, VARMAPhigran=NULL,
                        VARMAlThetaseq=NULL, VARMAThetagran=NULL, VARMAalpha=0,
                        VARMAselection = c("none", "cv", "bic", "aic", "hq"),
                        h=1, cvcut=0.9, eps=10^-3,
                        check_std = TRUE){

  ######################################
  #### Checking and Preparing Input ####
  ######################################
  VARselection <- match.arg(VARselection)
  VARMAselection <- match.arg(VARMAselection)


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
  if((!is.vector(VARMAlPhiseq) & !is.null(VARMAlPhiseq))){
    stop("The regularization parameter grid VARMAlPhiseq needs to be a vector of length >= 1 or NULL otherwise")
  }
  if (length(VARMAlPhiseq) == 1){
    if (length(VARMAlThetaseq) == 1 & VARMAselection != 'none'){
      stop("If both penalty sequences are of length one, then only selection='none' is supported")
    }
  }
  if(any((VARMAPhigran<=0)==T)){
    stop("The granularity parameters need to be a strictly positive integer")
  }
  if((!is.vector(VARMAlThetaseq) & !is.null(VARMAlThetaseq))){
    stop("The regularization parameter VARMAlThetaseq needs to be a vector of length >=1 or NULL otherwise")
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

  ###################################
  #### END Check and Preparation ####
  ###################################


  # Initialization
  HVARcvFIT  <- HVARFIT <-  NA
  HVARXcvFIT <- HVARXFIT <- NA
  VARPhi <- VARphi0 <- NA
  Phi <- Theta <- phi0 <- resid <- NA

  #######################
  #### START PHASE I ####
  #######################
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
  if (is.null(U)){
    #HVAR Estimation

    HVARFIT <- sparseVAR(Y = Y, p = VARp, VARlseq = VARlseq,
                      VARgran = c(VARgran1, VARgran2), selection = VARselection,
                      cvcut = cvcut, eps = eps, VARpen = VARpen,
                      check_std = FALSE)
    VARPhi <- HVARFIT$Phihat
    VARphi0 <- HVARFIT$phi0hat
    HVARresid <- residuals(HVARFIT)
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

  #####################
  #### END PHASE I ####
  #####################


  ########################
  #### START PHASE II ####
  ########################
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

  if (VARMAselection %in% c("bic", "aic", "hq")) warning("IC selection is not recommended for VARMA models.")
  HVARXFIT <- sparseVARX(Y = Y, X = U, p = VARMAp, s = VARMAq,
                         VARXpen = VARMApen, VARXlPhiseq = VARMAlPhiseq,
                         VARXPhigran = c(VARMAPhigran1, VARMAPhigran2),
                         VARXlBseq = VARMAlThetaseq,
                         VARXBgran = c(VARMAThetagran1, VARMAThetagran2),
                         VARXalpha = VARMAalpha, h = h, cvcut = cvcut,
                         eps = eps, selection = VARMAselection,
                         check_std = FALSE)

  Phi <- if (is.null(dim(HVARXFIT$Phihat))) matrix(HVARXFIT$Phihat, nrow = 1) else HVARXFIT$Phihat
  phi0 <- HVARXFIT$phi0hat
  Theta <- HVARXFIT$Bhat
  resid <- residuals(HVARXFIT)
  U <- matrix(U, ncol=ncol(U), nrow=nrow(U))
  k <- ncol(Y)

  ######################
  #### END PHASE II ####
  ######################

  # Output
  out<-list("k"=k, "Y"=Y, "VARp"=VARp, "VARPhihat"=VARPhi, "VARphi0hat"=VARphi0, "U"=U,
            "VARMAp"=VARMAp, "VARMAq"=VARMAq,
            "Phihat"=Phi, "Thetahat"=Theta, "phi0hat"=phi0,
            "series_names"=series_names,
            "PhaseI_lambdas"=HVARFIT$lambdas, "PhaseI_MSFEcv"=HVARFIT$MSFEcv,
            "PhaseI_lambda_SEopt"=HVARFIT$lambda_SEopt,"PhaseI_lambda_opt"=HVARFIT$lambda_opt,
            "PhaseII_lambdaPhi"=HVARXFIT$lambdaPhi, "PhaseII_lambdaTheta"=HVARXFIT$lambdaB,
            "PhaseII_lambdaPhi_opt"=HVARXFIT$lambdaPhi_opt, "PhaseII_lambdaPhi_SEopt"=HVARXFIT$lambdaPhi_SEopt,
            "PhaseII_lambdaTheta_opt"=HVARXFIT$lambdaB_opt, "PhaseII_lambdaTheta_SEopt"=HVARXFIT$lambdaB_SEopt,
            "PhaseII_MSFEcv"=HVARXFIT$MSFEcv, "h"=h,
            "VARselection" = VARselection, "VARMAselection" = VARMAselection)
  class(out) <- "bigtime.VARMA"
  out
}






