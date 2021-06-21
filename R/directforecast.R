
#' Function to obtain h-step ahead direct forecast based on estimated VAR, VARX or VARMA model
#' @param fit Fitted sparse VAR, VARX or VARMA model.
#' @param h Desired forecast horizon. Default is h=1.
#' @export
#' @return Vector of length k containing the h-step ahead forecasts for the k time series.
#' @examples
#' data(var.example)
#' VARfit <- sparseVAR(scale(Y.var), selection = "cv") # sparse VAR
#' VARforecast <- directforecast(fit=VARfit, h=1)
directforecast <- function(fit, h=1){

  model <- "none"
  if ("bigtime.VAR" %in% class(fit)) model <- "VAR"
  if ("bigtime.VARX" %in% class(fit)) model <- "VARX"
  if ("bigtime.VARMA" %in% class(fit)) model <- "VARMA"

  if(h<=0){
    stop("Forecast horizon h must be a strictly positive integer.")
  }

  if(fit$h!=h){
    stop("Specify the same forecast horizon h as argument in the sparseVAR, sparseVARX, or sparseVARMA function and the directforecast function")
  }

  if(!is.element(model, c("VAR", "VARX", "VARMA"))){
    stop("The model needs to be either VAR, VARX or VARMA")
  }

  # Preliminaries
  k <- ncol(fit$Y)

  if(model=="VAR"){

    # We need to catch the situation in which no selection was done
    if (fit$selection == "none"){
      fcst <- array(dim = c(1, ncol(fit$Y), length(fit$lambdas)))
      for (i in 1:length(fit$lambdas)){
        mod_tmp <- fit
        mod_tmp$selection <- "tmp"
        mod_tmp$Phihat <- fit$Phihat[, , i]
        mod_tmp$phi0hat <- fit$phi0hat[, , i]
        fcst[, , i] <- directforecast(mod_tmp, h = h)
      }
      return(fcst)
    }


    Y <- fit$Y
    p <- fit$p
    Phi <- fit$Phihat
    phi0 <- fit$phi0hat

    if(is.null(Phi) | is.null(phi0)){
      stop("Please provide a fitted VAR model")
    }

    VARFIT <- HVARmodel(Y=Y, p=p, h=h)

    if(k==1){
      VARFIT$fullY <- as.matrix(VARFIT$fullY)
      Zpred <- c(c(t(VARFIT$fullY[nrow(VARFIT$fullY):(nrow(VARFIT$fullY)-p+1),1:ncol(VARFIT$fullY)])))
      Phi <- matrix(Phi, nrow=1)
    }else{
      Zpred <- c(c(t(VARFIT$fullY[nrow(VARFIT$fullY):(nrow(VARFIT$fullY)-p+1),1:ncol(VARFIT$fullY)])))
    }
    Ypred <- matrix(Zpred, nrow=VARFIT$p, byrow=T)

    Yhat <- phi0
    for(i.p in 1:p){

      if(k==1){
        Yhat <- Yhat + matrix(Phi[, ((i.p-1)*k+1):(i.p*k)  ], nrow=k)%*%Ypred[i.p,]
      }else{
        Yhat <- Yhat + Phi[, ((i.p-1)*k+1):(i.p*k)  ]%*%Ypred[i.p,]
      }

    }
  }

  if(model=="VARX"| model=="VARMA"){


    if(model=="VARX"){
      m <- ncol(fit$X)
      Y <- fit$Y
      U <- fit$X
      p <- fit$p
      q <- fit$s
      Phi <- fit$Phihat
      Theta <- fit$Bhat
      phi0 <- fit$phi0hat

      if (fit$selection == "none"){
        dir_fcst <- array(dim = c(1, ncol(Y), dim(Phi)[3]))
        for (i in 1:dim(Phi)[3]){
          mod_tmp <- fit
          mod_tmp$Phihat <- fit$Phihat[,,i]
          mod_tmp$phi0hat <- fit$phi0hat[,,i]
          mod_tmp$Bhat <- fit$Bhat[,,i]
          mod_tmp$selection = "tmp"
          dir_fcst[,,i] <- directforecast(mod_tmp, h = h)
        }
        return(dir_fcst)
      }

      if(is.null(Phi) | is.null(Theta) | is.null(phi0)){
        stop("Please provide a fitted VARX model")
      }
    }

    if(model=="VARMA"){
      m <- ncol(fit$U)
      Y <- fit$Y
      U <- fit$U
      p <- fit$VARMAp
      q <- fit$VARMAq
      Phi <- fit$Phihat
      Theta <- fit$Thetahat
      phi0 <- fit$phi0hat

      if (fit$VARMAselection == "none"){
        # If no selection was done, then fitted values are tensors
        dir_fcst <- array(dim = c(1, ncol(Y), dim(Phi)[3]))
        for (i in 1:length(fit$PhaseII_lambdaPhi)){
          mod_tmp <- fit
          mod_tmp$Phihat <- fit$Phihat[, , i]
          mod_tmp$Thetahat <- fit$Thetahat[, , i]
          mod_tmp$phi0hat <- fit$phi0hat[, , i]
          mod_tmp$VARMAselection <- "tmp"
          dir_fcst[, , i] <- directforecast(mod_tmp, h = h)
        }
        return(dir_fcst)
      }

      if(is.null(Phi) | is.null(Theta) | is.null(phi0)){
        stop("Please provide a fitted VARMA model")
      }

    }

    VARXFIT <- HVARXmodelFORECAST(Y=Y, X=U, p=p, s=q, h=h)

    if(k==1){
      VARXFIT$fullY <- as.matrix(VARXFIT$fullY)
      Phi <- matrix(Phi, nrow=1)
      Theta <- matrix(Theta, nrow=1)
    }
    Zpred <- c(c(t(VARXFIT$fullY[nrow(VARXFIT$fullY):(nrow(VARXFIT$fullY)-p+1), 1:ncol(VARXFIT$fullY)])))
    Ypred <- matrix(Zpred, nrow=VARXFIT$p, byrow=T)

    if(m==1){
      VARXFIT$fullXRESP <- as.matrix(VARXFIT$fullXRESP)
    }
    Xpred <- c(c(t(VARXFIT$fullXRESP[nrow(VARXFIT$fullXRESP):(nrow(VARXFIT$fullXRESP)-q+1),1:ncol(VARXFIT$fullXRESP)])))
    Upred <- matrix(Xpred,nrow=VARXFIT$s,byrow=T)

    Yhat <- phi0
    for(i.p in 1:p){

      if(k==1){
        Yhat <- Yhat + matrix(Phi[, ((i.p-1)*k+1):(i.p*k)  ], nrow=k)%*%Ypred[i.p,]
      }else{
        Yhat <- Yhat + Phi[, ((i.p-1)*k+1):(i.p*k)  ]%*%Ypred[i.p,]
      }

    }
    for(i.q in 1:q){

      if(m==1){
        Yhat <- Yhat + matrix(Theta[, ((i.q-1)*m+1):(i.q*m)  ], nrow=k)%*%Upred[i.q,]
      }else{
        Yhat <- Yhat + Theta[, ((i.q-1)*m+1):(i.q*m)  ]%*%Upred[i.q,]
      }

    }

  }


  return("Yhat"=c(Yhat))
}

HVARXmodelFORECAST<-function(Y, X, p, s, h=1){

  # Preliminaries
  k <- ncol(Y)
  kX <-ncol(X)

  # Response and predictor matrices
  m <- max(s, p)
  DATAY <- embed(Y, dimension=m+h)
  fullY <- DATAY[, 1:k]

  fullZNEW <- as.matrix(as.matrix(DATAY[,-c(1:k)])[,(1:((p+h-1)*k))])
  fullZ <- t(fullZNEW[,(ncol(fullZNEW)-k*p+1):ncol(fullZNEW)])

  DATAX <- embed(X,dimension=m+h)
  fullXRESP <- DATAX[,1:kX]
  fullXNEW <- as.matrix(as.matrix(DATAX[,-c(1:kX)])[,(1:((s+h-1)*kX))])
  fullX <- t(fullXNEW[,(ncol(fullXNEW)-kX*s+1):ncol(fullXNEW)])

  out <- list("fullY"=fullY, "fullX"=fullX, "fullZ"=fullZ, "k"=k, "kX"=kX, "p"=p, "s"=s,
              "fullXRESP"=fullXRESP)
}
