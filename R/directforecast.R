
#' Function to obtain h-step ahead direct forecast based on estimated VAR, VARX or VARMA model
#' @param fit Fitted var, varx or varma model.
#' @param model Type of model that was estimated: var, varx or varma.
#' @param h Desired forecast horizon in time-series cross-validation procedure.
#' @export
#' @return A list with the following components
#' \item{Yhat}{h-step ahead forecast for the k time series}
#' @examples
#' data(Y)
#' varfit <- sparsevar(Y) # sparse VAR
#' forecasts <- directforecast(fit=varfit, model="var", h=1)
directforecast <- function(fit, model, h=1){
  # Preliminaries
  k <- ncol(fit$Y)

  if(model=="var"){
    Y <- fit$Y
    p <- fit$p
    Phi <- fit$Phihat
    phi0 <- fit$phi0hat

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

  if(model=="varx"| model=="varma"){


    if(model=="varx"){
      m <- ncol(fit$X)
      Y <- fit$Y
      U <- fit$X
      p <- fit$p
      q <- fit$s
      Phi <- fit$Phihat
      Theta <- fit$Bhat
      phi0 <- fit$phi0hat
    }

    if(model=="varma"){
      m <- ncol(fit$U)
      Y <- fit$Y
      U <- fit$U
      p <- fit$VARMAp
      q <- fit$VARMAq
      Phi <- fit$Phihat
      Theta <- fit$Thetahat
      phi0 <- fit$phi0hat
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


  return("Yhat"=Yhat)
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
