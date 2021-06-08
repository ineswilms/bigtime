#' Creates Lagmatrix of Estimated Coefficients
#' @param fit Fitted VAR, VARX or VARMA model.
#' @param returnplot TRUE or FALSE: return plot of lag matrix or not.
#' @export
#' @return A list with estimated lag matrix of the VAR model, or lag matrices of the VARX or VARMA model. The rows contain the responses, the columns contain the predictors.
#' @examples
#' data(Y)
#' data(X)
#' VARXfit <- sparseVARX(Y=scale(Y), X=scale(X), selection = "cv") # sparse VARX
#' Lhats <- lagmatrix(fit=VARXfit)
lagmatrix <- function(fit, returnplot=F){

  model <- "none"
  if ("bigtime.VAR" %in% class(fit)) model <- "VAR"
  if ("bigtime.VARX" %in% class(fit)) model <- "VARX"
  if ("bigtime.VARMA" %in% class(fit)) model <- "VARMA"

  if(!is.element(model, c("VAR", "VARX", "VARMA"))){
    stop("The model needs to be either VAR, VARX or VARMA")
  }

  if(model=="VAR"){

    if (fit$selection == "none") stop("Model did not use any selection procedure. Please use a selection procedure in sparseVAR or call ic_selection on model first.")

    k <- fit$k
    p <- fit$p
    coef <- fit$Phihat

    if(is.null(coef)){
      stop("Please provide a fitted VAR model")
    }
    if(k==1){
      coef <- matrix(coef, nrow=1)
    }

    Lhat <- matrix(NA, ncol=k, nrow=k)

    for(i.row in 1:k){ # responses
      for(j.col in 1:k){ # predictors
        if( length( which(coef[i.row, seq(from=j.col, length=p, by=k)]!=0) )==0 ){
          Lhat[i.row, j.col] <- 0
        }else{
          Lhat[i.row, j.col] <- max( which(coef[i.row, seq(from=j.col, length=p, by=k)]!=0) )
        }
      }
    }



    if(returnplot){

      if(k==1){
        stop("Plot only supported for multivariate time series case (k>1).")
      }
      # Check time series names
      if (is.null(fit$series_names))
        colnames(Lhat) <- paste0("Pred.", 1:ncol(Lhat))
      else
        colnames(Lhat) <- fit$series_names
      if (is.null(fit$series_names))
        rownames(Lhat) <- paste0("Resp.", 1:nrow(Lhat))
      else
        rownames(Lhat) <- fit$series_names
      par(mfrow=c(1,1))
      plotlagmat(Lhat)
    }

    Lhat <- list("LPhi"=Lhat)
  }

  if(model=="VARX"| model=="VARMA"){

    if(model=="VARX"){
      k <- fit$k
      m <- fit$m
      p <- fit$p
      s <- fit$s
      coef1 <- fit$Phihat
      coef2 <- fit$Bhat

      if(fit$selection == "none") stop("No selection procedure was used.")
      if(is.null(coef1)|is.null(coef2)){
        stop("Please provide a fitted VARX model")
      }

    }

    if(model=="VARMA"){
      k <- fit$k
      m <- fit$k
      p <- fit$VARMAp
      s <- fit$VARMAq
      coef1 <- fit$Phihat
      coef2 <- fit$Thetahat

      if (fit$VARMAselection == "none") stop("No selection procedure was used.")
      if(is.null(coef1)|is.null(coef2)){
        stop("Please provide a fitted VARMA model")
      }
    }


    if(k==1){
      coef1 <- matrix(coef1, nrow=1)
      coef2 <- matrix(coef2, nrow=1)
    }

    Lhat1 <- matrix(NA, ncol=k, nrow=k)
    Lhat2 <- matrix(NA, ncol=m, nrow=k)

    for(i.row in 1:k){ # responses
      for(j.col in 1:k){ # predictors
        if( length( which(coef1[i.row, seq(from=j.col, length=p, by=k)]!=0) )==0 ){
          Lhat1[i.row, j.col] <- 0
        }else{
          Lhat1[i.row, j.col] <- max( which(coef1[i.row, seq(from=j.col, length=p, by=k)]!=0) )
        }
      }
    }

    for(i.row in 1:k){ # responses
      for(j.col in 1:m){ # predictors
        if( length( which(coef2[i.row, seq(from=j.col, length=s, by=m)]!=0) )==0 ){
          Lhat2[i.row, j.col] <- 0
        }else{
          Lhat2[i.row, j.col] <- max( which(coef2[i.row, seq(from=j.col, length=s, by=m)]!=0) )
        }
      }
    }

    if(model=="VARX"){
      Lhat <- list("LPhi"=Lhat1, "LB"= Lhat2)
    }

    if(model=="VARMA"){
      Lhat <- list("LPhi"=Lhat1, "LTheta"=Lhat2)
    }

    if(returnplot){

      if(k==1){
        warning("Plot only supported for multivariate time series case (k>1).")
        return(Lhat)
      }

      if(model=="VARX"){
        if (is.null(fit$endogenous_series_names))
          endog_names <- paste0("Endog.", 1:ncol(Lhat1))
        else
          endog_names <- fit$endogenous_series_names
        colnames(Lhat1) <- endog_names
        rownames(Lhat1) <- endog_names
        rownames(Lhat2) <- endog_names
        if (is.null(fit$exogenous_series_names))
          colnames(Lhat2) <- paste0("Exog.", 1:ncol(Lhat2))
        else
          colnames(Lhat2) <- fit$exogenous_series_names
      }

      if(model=="VARMA"){
        if (is.null(fit$series_names)) {
          colnames(Lhat1) <- paste0("AR.", 1:ncol(Lhat1))
          colnames(Lhat2) <- paste0("MA.", 1:ncol(Lhat2))
          rownames(Lhat1) <- paste0("Resp.", 1:nrow(Lhat1))
          rownames(Lhat2) <- paste0("Resp.", 1:nrow(Lhat2))
        }
        else {
          colnames(Lhat1) <- paste0("AR.", fit$series_names)
          colnames(Lhat2) <- paste0("MA.", fit$series_names)
          rownames(Lhat1) <- fit$series_names
          rownames(Lhat2) <- fit$series_names
        }
      }


      plotlagmat(Lhat1)
      plotlagmat(Lhat2)
    }
  }


  return(Lhat)
}


plotlagmat <- function(lagmat){
  Response <- Predictor <- Lags <- NULL # Needed because otherwise cran complaines about lazy evaluation


  lmat <- tidyr::as_tibble(lagmat)
  lmat <- lmat %>%
    dplyr::mutate(Response = rownames(lagmat)) %>%
    tidyr::pivot_longer(-Response, names_to = "Predictor", values_to = "Lags") %>%
    mutate(Response = factor(Response, levels = sort(unique(Response), decreasing = TRUE)),
           Predictor = factor(Predictor, levels = sort(unique(Predictor), decreasing = FALSE)))

  p <- ggplot2::ggplot(data = lmat,
                  ggplot2::aes(x = Predictor, y = Response, fill = Lags)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(x = Predictor, y = Response, label = Lags),
                       color = "black", size = 5) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_gradient(low = "#ffffff", high = "#147aff") +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x.top = ggplot2::element_text(margin = ggplot2::margin(b = 20),
                                               size = 15),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 20),
                                           size = 15),
      axis.text.x = ggplot2::element_text(face = "bold",
                                          size = 10, angle = 90),
      axis.text.y = ggplot2::element_text(face = "bold", size = 10)
    )

  plot(p)
}
