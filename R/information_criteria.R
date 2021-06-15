
#' Selects the optimal penalty parameter using information criteria
#'
#' @param mod Model estimated Model estimated using \code{\link{sparseVAR}},
#' \code{\link{sparseVARX}}, or \code{\link{sparseVARMA}}
#' @param ic Which information criteria should be used. Must be one of
#' \code{"bic"}, \code{"aic"} or \code{"hq"}
#' @param verbose If true, some useful information will be printed during the process
#' @export
#' @return Returns a model that uses the optimal penalty
ic_selection <- function(mod, ic = c("bic", "aic", "hq"), verbose = FALSE){
  ics <- get_ic_vals(mod, verbose = verbose)
  ic <- match.arg(ic)
  ic <- toupper(ic) # This just makes it more convenient. typing lower case is faster
  selected <- ics$mins[ic]
  obtain_selected_model(mod, selected)
}

# Creates a model out of another model and a selected penalisation
obtain_selected_model <- function(mod, selected){
  mod_new <- mod
  if ("bigtime.VAR" %in% class(mod)){
    Phihat <- mod$Phihat[, , selected]
    phi0hat <- mod$phi0hat[, , selected]
    mod_new$Phihat <- Phihat
    mod_new$phi0hat <- phi0hat
    mod_new$lambda_opt <- mod$lambdas[selected]
  }else if ("bigtime.VARX" %in% class(mod)){
    mod_new$Phihat <- mod$Phihat[, , selected]
    mod_new$phi0hat <- mod$phi0hat[, , selected]
    mod_new$Bhat <- mod$Bhat[, , selected]
    mod_new$lambdaPhi_opt <- mod$lambdaPhi[selected]
    mod_new$lambdaB_opt <- mod$lambdaB[selected]
  }else{
    stop("Unsopported model")
  }
  return(mod_new)
}

#' Calculates the Information Criteria for a VAR, VARX, VARMA model
#'
#' The number of non-zero coefficients are taken as the degrees of freedom.
#' This is not valid in all cases.
#'
#' @param mod Model estimated Model estimated using \code{\link{sparseVAR}},
#' \code{\link{sparseVARX}}, or \code{\link{sparseVARMA}}
#' @param verbose Should information about the optimal selection be printed?
#' @export
#'
#' @examples
#' dat <- simVAR(200, 2, 5, decay = 0.01)
#' mod <- sparseVAR(scale(dat$Y))
#' ics <- get_ic_vals(mod)
get_ic_vals <- function(mod, verbose = TRUE){
  UseMethod("get_ic_vals", mod)
}


#' Calculates the Information Criteria for a model estimated using
#' \code{\link{sparseVARX}}
#'
#' The number of non-zero coefficients in both the \code{Phihat}
#' and \code{Bhat} matrix are taken as the degrees of freedom.
#'
#' @param mod Model estimated using \code{\link{sparseVARX}}
#' @param verbose Should information about the optimal selection be printed?
#' @export
#' @return Returns a list containing
#' \item{ics}{Values of the ICs for all lambdas}
#' \item{mins}{Which IC lead to the minimum (the row number)}
#' \item{selected_lamPhi}{Which lambda Phi were selected}
#' \item{selected_lamB}{Which lambda B were selected}
get_ic_vals.bigtime.VARX <- function(mod, verbose = TRUE){
  if (!("bigtime.VARX" %in% class(mod))) stop("Currently only implemented for VARX models")
  tt <- nrow(mod$Y)
  if (mod$p*mod$k + mod$s*ncol(mod$X) > tt) warning("IC selection is not recommended.")
  if (is.matrix(mod$Phihat)) stop("Model was estimated using only one penalty parameter.")

  n3 <- dim(mod$Phihat)[3]
  ics <- matrix(ncol = 3, nrow = n3)
  mod_tmp <- mod
  mod_tmp$selection <- "none"
  res <- residuals.bigtime.VARX(mod_tmp)
  for (i in 1:n3){
    df <- sum(mod$Phihat != 0) + sum(mod$Bhat != 0)
    r <- matrix(res[, , i], nrow = dim(res)[1], ncol = dim(res)[2])
    ics[i, ] <- .get_ic_vals(r, tt, df)
  }
  colnames(ics) <-  c("AIC", "BIC", "HQ")
  lambdaPhi <- mod$lambdaPhi
  lambdaB <- mod$lambdaB
  ics_df <- as.data.frame(cbind(lambdaPhi, lambdaB, ics))
  mins <- apply(ics, 2, which.min)
  selected_lamPhi <- do.call(c, lapply(mins, function(x) lambdaPhi[[x]]))
  selected_lamB <- do.call(c, lapply(mins, function(x) lambdaB[[x]]))
  names(selected_lamPhi) <- names(selected_lamB) <- names(mins)

  if (verbose){
    cat("\n\n#### Selected the following ####\n\n")
    tmp <- rbind(selected_lamPhi, selected_lamB)
    colnames(tmp) <- colnames(ics)
    rownames(tmp) <- c("lambda Phi:", "lambda B:")
    print(tmp)
    # cat("lambda Phi: \t", selected_lamPhi, "\n")
    # cat("lambda B: \t", selected_lamB, "\n")

    selected_aic <- rep("   ", nrow(ics))
    selected_aic[mins[1]] <- "==>"
    selected_bic <- rep("   ", nrow(ics))
    selected_bic[mins[2]] <- "==>"
    selected_hq <- rep("   ", nrow(ics))
    selected_hq[mins[3]] <- "==>"

    decimals <- 4
    aic_text <- paste(selected_aic, round(ics[, 1], decimals))
    bic_text <- paste(selected_bic, round(ics[, 2], decimals))
    hq_text <- paste(selected_hq, round(ics[, 3], decimals))
    lamPhi_text <- as.character(round(lambdaPhi, decimals))
    lamB_text <- as.character(round(lambdaB, decimals))
    text <- cbind(lamPhi_text, lamB_text, aic_text, bic_text, hq_text)
    colnames(text) <- c("lambda Phi", "lambda B", "AIC", "BIC", "HQ")

    cat("\n\n#### Details ####\n\n")
    print(as.data.frame(text))
  }

  list(
    ics = ics,
    mins = mins,
    selected_lamPhi = selected_lamPhi,
    selected_lamB = selected_lamB
  )

}

#' Calculates the Information Criteria for a model estimated using
#' \code{\link{sparseVAR}}
#'
#' The number of non-zero coefficients are taken as the degrees of freedom.
#'
#' @param mod Model estimated using \code{\link{sparseVAR}}
#' @param verbose Should information about the optimal selection be printed?
#' @export
#' @return Returns a list containing
#' \item{ics}{Values of the ICs for all lambdas}
#' \item{mins}{Which IC lead to the minimum (the row number)}
#' \item{selected_lambdas}{Which lambdas were selected}
#'
#' @examples
#' dat <- simVAR(200, 2, 5, decay = 0.01)
#' mod <- sparseVAR(scale(dat$Y))
#' ics <- get_ic_vals(mod)
get_ic_vals.bigtime.VAR <- function(mod, verbose = TRUE){
  if (!("bigtime.VAR" %in% class(mod))) stop("Currently only implemented for VAR models")
  tt <- nrow(mod$Y)
  if (mod$p * mod$k >= tt) warning("IC selection is only recommended for cases in which p < T")
  if (is.matrix(mod$Phihat)) stop("Model was estimated using only one penalty parameter.")

  n3 <- dim(mod$Phihat)[3]
  ics <- matrix(ncol = 3, nrow = n3)
  for (i in 1:n3){
    mod_tmp <- mod
    mod_tmp$Phihat <- mod$Phihat[, , i]
    mod_tmp$phi0hat <- mod$phi0hat[, , i]
    mod_tmp$selection <- "tmp"

    res <- residuals(mod_tmp)
    df <- sum(mod$Phihat[, , i] != 0)
    ics[i, ] <- .get_ic_vals(res, tt, df)
  }
  colnames(ics) <- c("AIC", "BIC", "HQ")
  lambda <- mod$lambdas
  ics_df <- as.data.frame(cbind(lambda, ics))
  mins <- apply(ics, 2, which.min)
  selected_lam <- do.call(c, lapply(mins, function(x) lambda[[x]]))
  names(selected_lam) <- names(mins)

  if (verbose){
    cat("\n\n#### Selected the following lambda ####\n\n")
    print(selected_lam)

    selected_aic <- rep("   ", nrow(ics))
    selected_aic[mins[1]] <- "==>"
    selected_bic <- rep("   ", nrow(ics))
    selected_bic[mins[2]] <- "==>"
    selected_hq <- rep("   ", nrow(ics))
    selected_hq[mins[3]] <- "==>"

    decimals <- 4
    aic_text <- paste(selected_aic, round(ics[, 1], decimals))
    bic_text <- paste(selected_bic, round(ics[, 2], decimals))
    hq_text <- paste(selected_hq, round(ics[, 3], decimals))
    lam_text <- as.character(round(lambda, decimals))
    text <- cbind(lam_text, aic_text, bic_text, hq_text)
    colnames(text) <- c("lambda", "AIC", "BIC", "HQ")

    cat("\n\n#### Details ####\n\n")
    print(as.data.frame(text))
  }

  list(
    ics = ics_df,
    mins = mins,
    selected_lambda = selected_lam
  )
}

#' Calculates ICs
#'
#' Not meant to be called by the user
#'
#' @param res Matrix of residuals
#' @param tt Number of time periods
#' @param df Degrees of freedom
#' @keywords internal
.get_ic_vals <- function(res, tt, df){

  Omega <- cov(res)
  log_det_omega <- log(det(Omega))

  aic <- log_det_omega + 2/tt * df
  bic <- log_det_omega + log(tt)/tt * df
  hq <- log_det_omega + 2*log(log(tt))/tt * df

  out <- cbind(aic, bic, hq)

  out
}


