#' Creates lagmatrix of estimated coefficients
#' @param fit Fitted var, varx or varma model.
#' @param model Type of model that was estimated: var, varx or varma.
#' @param returnplot True or False: return plot of Lhat or not.
#' @export
#' @return A list with the following components
#' \item{Lhat}{Estimated lag matrix in the VAR model, or lag matrices in the VARX or VARMA model. The rows contain the responses, the columns contain the predictors.}
#' @examples
#' data(Y)
#' data(X)
#' varxfit <- sparsevarx(Y=Y, X=X) # sparse varx
#' lhats <- lagmatrix(fit=varxfit, model="varx")
lagmatrix <- function(fit, model, returnplot=F){

  if(!is.element(model, c("var", "varx", "varma"))){
    stop("The model needs to be either var, varx or varma")
  }

  if(model=="var"){
    k <- fit$k
    p <- fit$p
    coef <- fit$Phihat

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
      colnames(Lhat) <- paste0("Pred.", 1:ncol(Lhat))
      rownames(Lhat) <- paste0("Resp.", 1:nrow(Lhat))
      par(mfrow=c(1,1))
      plotlaghat(datamatrix=Lhat, cl.lim=c(0, p),
                 title="Lhat", mar=c(0.5, 0.1, 2, 0.1))
    }

    Lhat <- list("LPhi"=Lhat)
  }

  if(model=="varx"| model=="varma"){

    if(model=="varx"){
      k <- fit$k
      m <- fit$m
      p <- fit$p
      s <- fit$s
      coef1 <- fit$Phihat
      coef2 <- fit$Bhat
    }

    if(model=="varma"){
      k <- fit$k
      m <- fit$k
      p <- fit$VARMAp
      s <- fit$VARMAq
      coef1 <- fit$Phihat
      coef2 <- fit$Thetahat
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

    if(model=="varx"){
      Lhat <- list("LPhi"=Lhat1, "LB"= Lhat2)
    }

    if(model=="varma"){
      Lhat <- list("LPhi"=Lhat1, "LTheta"=Lhat2)
    }

    if(returnplot){

      if(k==1){
        warning("Plot only supported for multivariate time series case (k>1).")
        return(Lhat)
      }

      if(model=="varx"){
        colnames(Lhat1) <- paste0("Endog.", 1:ncol(Lhat1))
        colnames(Lhat2) <- paste0("Exog.", 1:ncol(Lhat2))
      }

      if(model=="varma"){
        colnames(Lhat1) <- paste0("AR.", 1:ncol(Lhat1))
        colnames(Lhat2) <- paste0("MA.", 1:ncol(Lhat2))
      }

      rownames(Lhat1) <- paste0("Resp.", 1:nrow(Lhat1))
      rownames(Lhat2) <- paste0("Resp.", 1:nrow(Lhat2))

      plotlaghat(datamatrix=Lhat1, cl.lim=c(0, p),
                 title="Lhat", mar=c(0.5, 0.1, 2, 0.1))

      plotlaghat(datamatrix=Lhat2, cl.lim=c(0, s),
                 title="Lhat", mar=c(0.5, 0.1, 2, 0.1))
    }
  }


  return(Lhat)
}


plotlaghat <- function (datamatrix, title = "", mar = c(0, 0, 0, 0), cl.lim, addlegend=T,
                            norowname=F, nocolname=F)
{ # This function is based on the corrplot function from the corrplot Rpackage

  corr <- datamatrix

  # Settings from corrplot function
  type <- "full"
  outline <- T
  diag <- T
  addCoefasPercent <- FALSE
  tl.offset <- 0.4
  tl.srt <- 90
  tl.cex <- 1.2
  tl.col <- "black"
  cl.ratio <- 0.1
  cl.length <- NULL
  cl.align.text <- "c"
  cl.offset <- 0.5
  cl.cex <- 0.8
  number.cex <- 1
  number.font <- 2
  number.digits <- NULL
  bg <- "white"
  addgrid.col <- "grey"
  addCoef.col <- "black"
  rect.col <-  "black"


  intercept <- -cl.lim[1]
  zoom <- 1/(diff(cl.lim))
  corr <- (intercept + corr) * zoom

  cl.lim2 <- (intercept + cl.lim) * zoom
  int <- intercept * zoom
  col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D",
                            "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                            "#4393C3", "#2166AC", "#053061"))(200)
  n <- nrow(corr)
  m <- ncol(corr)
  min.nm <- min(n, m)
  ord <- seq_len(min.nm)
  apply_mat_filter <- function(mat) {
    x <- matrix(1:n * m, n, m)
    switch(type, upper = mat[row(x) > col(x)] <- Inf, lower = mat[row(x) <
                                                                    col(x)] <- Inf)
    if (!diag) {
      diag(mat) <- Inf
    }
    return(mat)
  }
  getPos.Dat <- function(mat) {
    tmp <- apply_mat_filter(mat)
    Dat <- tmp[is.finite(tmp)]
    ind <- which(is.finite(tmp), arr.ind = TRUE)
    Pos <- ind
    Pos[, 1] <- ind[, 2]
    Pos[, 2] <- -ind[, 1] + 1 + n
    return(list(Pos, Dat))
  }
  Pos <- getPos.Dat(corr)[[1]]
  n2 <- max(Pos[, 2])
  n1 <- min(Pos[, 2])
  nn <- n2 - n1
  newrownames <- as.character(rownames(corr)[(n + 1 - n2):(n +
                                                             1 - n1)])
  m2 <- max(Pos[, 1])
  m1 <- min(Pos[, 1])
  mm <- max(1, m2 - m1)
  newcolnames <- as.character(colnames(corr)[m1:m2])
  DAT <- getPos.Dat(corr)[[2]]
  len.DAT <- length(DAT)
  assign.color <- function(dat = DAT, color = col) {
    newcorr <- (dat + 1)/2
    newcorr[newcorr <= 0] <- 0
    newcorr[newcorr >= 1] <- 1 - 1e-16
    color[floor(newcorr * length(color)) + 1]
  }
  col.fill <- assign.color()
  isFALSE <- function(x) identical(x, FALSE)
  isTRUE <- function(x) identical(x, TRUE)
  tl.pos <- switch(type, full = "lt", lower = "ld", upper = "td")
  cl.pos <- switch(type, full = "r", lower = "b", upper = "r")
  col.border <- "black"

  oldpar <- par(mar = mar, bg = "white")
  on.exit(par(oldpar), add = TRUE)

  plot.new()
  xlabwidth <- ylabwidth <- 0
  for (i in 1:50) {
    xlim <- c(m1 - 0.5 - xlabwidth, m2 + 0.5 + mm * cl.ratio *
                (cl.pos == "r"))
    ylim <- c(n1 - 0.5 - nn * cl.ratio * (cl.pos == "b"),
              n2 + 0.5 + ylabwidth)
    plot.window(xlim + c(-0.2, 0.2), ylim + c(-0.2, 0.2),
                asp = 1, xaxs = "i", yaxs = "i")
    x.tmp <- max(strwidth(newrownames, cex = tl.cex))
    y.tmp <- max(strwidth(newcolnames, cex = tl.cex))
    if (min(x.tmp - xlabwidth, y.tmp - ylabwidth) < 1e-04) {
      break
    }
    xlabwidth <- x.tmp
    ylabwidth <- y.tmp
  }
  if (tl.pos == "n" || tl.pos == "d") {
    xlabwidth <- ylabwidth <- 0
  }
  if (tl.pos == "td")
    ylabwidth <- 0
  if (tl.pos == "ld")
    xlabwidth <- 0
  laboffset <- strwidth("W", cex = tl.cex) * tl.offset
  xlim <- c(m1 - 0.5 - xlabwidth - laboffset, m2 + 0.5 +
              mm * cl.ratio * (cl.pos == "r")) + c(-0.35, 0.15)
  ylim <- c(n1 - 0.5 - nn * cl.ratio * (cl.pos == "b"),
            n2 + 0.5 + ylabwidth * abs(sin(tl.srt * pi/180)) +
              laboffset)
  +c(-0.15, 0.35)
  if (.Platform$OS.type == "windows") {
    grDevices::windows.options(width = 7, height = 7 *
                                 diff(ylim)/diff(xlim))
  }
  plot.window(xlim = xlim, ylim = ylim, asp = 1, xlab = "",
              ylab = "", xaxs = "i", yaxs = "i")

  laboffset <- strwidth("W", cex = tl.cex) * tl.offset
  symbols(Pos, add = TRUE, inches = FALSE, squares = rep(1,
                                                         len.DAT), bg = bg, fg = bg)


  number.digits <- switch(addCoefasPercent + 1, 2, 0)

  stopifnot(number.digits%%1 == 0)
  stopifnot(number.digits >= 0)

  NA_LABEL_MAX_CHARS <- 2


  symbols(Pos, add = TRUE, inches = FALSE, squares = rep(1,
                                                         len.DAT), bg = col.fill, fg = col.border)

  symbols(Pos, add = TRUE, inches = FALSE, bg = NA, squares = rep(1, len.DAT), fg = addgrid.col)

  colRange <- assign.color(dat = cl.lim2)
  ind1 <- which(col == colRange[1])
  ind2 <- which(col == colRange[2])
  colbar <- col[ind1:ind2]
  if (is.null(cl.length)) {
    cl.length <- ifelse(length(colbar) > 20, 11, length(colbar) +
                          1)
  }
  labels <- seq(cl.lim[1], cl.lim[2], length = cl.length)
  if (cl.pos == "r") {
    vertical <- TRUE
    xlim <- c(m2 + 0.5 + mm * 0.02, m2 + 0.5 + mm * cl.ratio)
    ylim <- c(n1 - 0.5, n2 + 0.5)
  }
  if (cl.pos == "b") {
    vertical <- FALSE
    xlim <- c(m1 - 0.5, m2 + 0.5)
    ylim <- c(n1 - 0.5 - nn * cl.ratio, n1 - 0.5 - nn *
                0.02)
  }

  if(addlegend){
    # Add legend
    colorlegend(colbar = colbar, labels = round(labels, 0),
                offset = cl.offset, ratio.colbar = 0.3, cex = cl.cex,
                xlim = xlim, ylim = ylim, vertical = vertical, align = cl.align.text)
  }


  pos.xlabel <- cbind(m1:m2, n2 + 0.5 + laboffset)
  pos.ylabel <- cbind(m1 - 0.5, n2:n1)

  if(norowname){
    newrownames <- rep("",length(newrownames))
  }

  if(nocolname){
    newcolnames <- rep("",length(newcolnames))
  }

  text(pos.xlabel[, 1], pos.xlabel[, 2], newcolnames,
       srt = tl.srt, adj = ifelse(tl.srt == 0, c(0.5,
                                                 0), c(0, 0)), col = tl.col, cex = tl.cex, offset = tl.offset)
  text(pos.ylabel[, 1], pos.ylabel[, 2], newrownames,
       col = tl.col, cex = tl.cex, pos = 2, offset = tl.offset)


  # title(title, ...)

  # Add numbers to each non-zero cell
  check <- c(datamatrix)
  check[check==0] <- c("")
  text(Pos[, 1], Pos[, 2], col = addCoef.col, labels = check,
       cex = number.cex, font = number.font)


}
