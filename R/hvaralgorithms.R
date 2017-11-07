# Elementwise HVAR
.HVARElemAlg <- function (beta, Y, Z, lambda, eps,p,MN,C1)
{

  k <- ncol(Y)
  if(MN)
  {
    C <- matrix(0,nrow=k,ncol=k*p)
    ## diag(C) <- rep(1,k)
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }

  k <- ncol(Y)
  betafin <- beta
  YOLD <- Y
  ZOLD <- Z
  YMean <- c(apply(Y, 2, mean))
  if(k==1 & p==1){
    ZMean <- matrix(mean(Z), ncol=1, nrow=1)
  }else{
    ZMean <- c(apply(Z, 1, mean))
  }


  Y <- Y - matrix(c(rep(1, nrow(Y))), ncol = 1) %*% matrix(c(apply(Y,
                                                                   2, mean)), nrow = 1)
  Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))

  tk <- 1/max(Mod(eigen(Z %*% t(Z))$values))


  if(k==1){
    betaini <- array(beta[, 2:length(beta[, , 1]), ],dim=c(k,k*p,length(lambda)))
    if(p==1){
      Z <- matrix(Z, nrow=1)
    }
  }else{
    betaini <- array(beta[, 2:ncol(beta[, , 1]), ],dim=c(k,k*p,length(lambda)))
  }



  betafin <- gamloopElem(betaini,Y,Z,lambda,eps,YMean,ZMean,as.matrix(betaini[,,1]),k,p)

  if(MN)
  {
    for(i in 1:(dim(betafin)[3]))
    {
      if(k==1){
        betafin[,2:length(betafin[,,i]),i] <- betafin[,2:length(betafin[,,i]),i]+C
      }else{
        betafin[,2:ncol(betafin[,,i]),i] <- betafin[,2:ncol(betafin[,,i]),i]+C
      }

    }

  }

  return(betafin)
}




# BigVARSupportFunctions

.LambdaGridE<- function (gran1, gran2, jj = jj, Y, Z, group,p,k,MN,alpha,C)
{


  if (group == "Basic"|group=="Tapered") {

    gamstart <- max(t(Y) %*% t(Z))

  }


  if (group == "HVARC"|group=="HVAROO"|group=="HVARELEM"|group=="HVARELEMSQRTLAS"|group=="HVARELEMnucl") {

    gmax <- c()

    for (i in 1:k) {

      gmax[i]  <- norm2(Z %*% Y[, i])

    }

    gamstart <- max(gmax)

  }

  if(group=="Tapered"){

    beta <- array(0,dim=c(k,k*p+1,10))

  }else{
    beta <- array(0,dim=c(k,k*p+1,1))
  }

  gamstart <- LGSearch(gamstart,Y,Z,beta,group,k,p,jj,MN,alpha,C)

  gamm <- exp(seq(from = log(gamstart), to = log(gamstart/gran1),
                  length = gran2))

  return(gamm)

}


LGSearch <- function(gstart,Y,Z,BOLD,group,k,p,gs,MN,alpha,C)
{

  tk <- 1/max(Mod(eigen(Z%*%t(Z))$values))
  lambdah <- gstart
  lambdal <- 0
  activeset <- list(rep(rep(list(0), length(gs))))



  while(max(abs(lambdah-lambdal))>.00001)
  {

    lambda <- (lambdah+lambdal)/2
    if(group=="Basic"){
      param <- .lassoVARFist(BOLD,Z,Y,lambda,1e-04,p,MN,C)[,2:(k*p+1),]
    }



    if(group=="HVARELEM")
    {
      BOLD <- .HVARElemAlg(BOLD,Y,Z,lambda,1e-4,p,MN,C)
      param <- BOLD[,2:(k*p+1),]
    }


    if(MN){

      diag(param[1:k,1:k]) <- ifelse(C==0,diag(param[1:k,1:k]),0)
      diag(BOLD[,2:(k*p+1),]) <- ifelse(C==0,diag(BOLD[,2:(k*p+1),]),0)

    }

    if(max(abs(param))==0)
    {
      lambdah <- lambda

    }else{
      lambdal <- lambda


    }


  }

  lambdah

}



.lfunction3 <- function(p,k)
{

  kk <- .lfunction2(p,k)

  oo <- list()

  pp <- list()

  for(i in 1:length(kk))
  {
    j <- 0
    oo[[i]] <- kk[[i]][(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
    pp[[i]] <- kk[[i]][-(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
    j <- j+1

  }

  ownoth <- c(oo,pp)
  return(ownoth)
}


.lfunction2 <- function(p,k)
{

  kk <- list()

  kk[[1]] <- 1:(k^2)

  if(p>1)
  {
    for(i in 2:p)
    {
      kk[[i]] <- 1:(k^2)+tail(kk[[i-1]],1)

    }
  }
  return(kk)
}


#Basic VAR Fista Implementation
.lassoVARFist <- function (B, Z, Y, gamm, eps,p,MN,C1)
{

  if(class(Y)!="matrix")
  {
    Y <- matrix(Y,ncol=1)
  }

  k <- ncol(Y)

  if(MN)
  {
    C <- matrix(0,nrow=k,ncol=k*p)
    diag(C) <- C1
    Y <- t(Y)
    Y <- Y-C%*%Z
    Y <- t(Y)
  }

  Y <- t(Y)
  YMean <- c(apply(Y, 1, mean))
  if(k==1 & p==1){
    ZMean <- matrix(mean(Z), ncol=1, nrow=1)
  }else{
    ZMean <- c(apply(Z, 1, mean))
  }
  Y <- Y - YMean %*% t(c(rep(1, ncol(Y))))


  if(k==1 & p==1){
    Z <- Z - c(apply(Z, 1, mean)) %*% t(c(rep(1, ncol(Z))))
  }else{
    Z <- Z - ZMean %*% t(c(rep(1, ncol(Z))))
  }


  Y <- t(Y)
  tk <- 1/max(Mod(eigen(Z%*%t(Z),only.values=TRUE)$values))

  BFOO1 <- matrix(B[, 2:dim(B)[2], 1],nrow=k,ncol=k*p)
  BFOO <- array(B[,2:ncol(as.matrix(B[,,1])),],dim=c(k,k*p,length(gamm)))



  if(k==1 & p==1){
      Z <- matrix(Z, nrow=1)
  }


  beta <- gamloopFista(BFOO, Y, Z, gamm, eps,
                       as.matrix(YMean), as.matrix(ZMean), BFOO1,k,p,tk,k,p)


  if(MN)
  {
    for(i in 1:(dim(beta)[3]))
      beta[,2:ncol(beta[,,i]),i] <- beta[,2:ncol(beta[,,i]),i]+C
  }
  return(beta)
}
