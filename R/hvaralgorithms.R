.lfunction3 <- function(p, k)
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


.lfunction2 <- function(p, k)
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


.LambdaGridE <- function (gran1, gran2, jj = jj, Y, Z, group, p, k, MN, alpha, C)
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

  if (group == "HVARELEM" | group =="Basci"){
    estim <- ifelse(group == "Basic", 1, 2)
    gamstart <- LGSearch_cpp(gamstart, Y, Z, beta, estim, k, p)
  }
  else{
    gamstart <- LGSearch(gamstart,Y,Z,beta,group,k,p,jj,MN,alpha,C)
  }

  gamm <- exp(seq(from = log(gamstart), to = log(gamstart/gran1), length = gran2))
  return(gamm)
}


LGSearch<- function(gstart, Y, Z, BOLD, group, k, p, gs, MN, alpha, C)
{

  # tk <- 1/max(Mod(eigen(Z%*%t(Z))$values))
  lambdah <- gstart
  lambdal <- 0
  activeset <- list(rep(rep(list(0), length(gs))))


  while(max(abs(lambdah-lambdal))>.00001)
  {

    lambda <- (lambdah+lambdal)/2


    if(group=="Basic"){
      param <- lassoVARFistcpp(BOLD, Y, Z, lambda, 1e-04, p)[,2:(k*p+1),]
    }



    if(group=="HVARELEM")
    {
      BOLD <- HVARElemAlgcpp(BOLD, Y, Z, lambda, 1e-4, p)
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
