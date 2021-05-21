
#' Selects the optimal penalty parameter using information criteria
#'
#' @param mod Model estimated using sparseVAR
#' @param ic Which information criteria should be used. Must be one of BIC, AIC, HQ
#' @param verbose If true, some useful information will be printed during the process
#' @export
#' @return returns a bigtime.VAR model that uses the optimal penalty
ic_selection <- function(mod, ic = c("bic", "aic", "hq"), verbose = FALSE){
  if (!("bigtime.VAR" %in% class(mod))) stop("Currently only implemented for VAR models")
  ics <- get_ic_vals(mod, verbose = verbose)
  ic <- match.arg(ic)
  ic <- toupper(ic) # This just makes it more convenient. typing lower case is faster
  selected <- ics$mins[ic]
  Phihat <- mod$Phihat[, , selected]
  phi0hat <- mod$phi0hat[, , selected]
  mod$Phihat <- Phihat
  mod$phi0hat <- phi0hat
  mod$lambda_opt <- mod$lambdas[selected]
  mod
}


#' Calculates the Information Criteria for a model estimated using sparseVAR
#'
#' The number of non-zero coefficients are taken as the degrees of freedom.
#' This is not valid in cases in which p>N
#'
#' @param mod Model estimated using sparseVAR
#' @param verbose Should information about the optimal selection be printed?
#' @export
#' @return Returns a list containing
#' \item{ics}{Values of the ICs for all lambdas}
#' \item{mins}{Which ic lead to the minimum (the row number)}
#' \item{selected_lambdas}{Which lambdas were selected}
#'
#' @examples
#' dat <- bigtime::simVAR(200, 2, 5, decay = 0.01)
#' mod <- bigtime::sparseVAR(scale(dat$Y), cv=FALSE)
#' ics <- get_ic_vals(mod)
get_ic_vals <- function(mod, verbose = TRUE){
  if (!("bigtime.VAR" %in% class(mod))) stop("Currently only implemented for VAR models")
  tt <- nrow(mod$Y)
  if (mod$p * mod$k >= tt) warning("IC selection is only recommended for cases in which p < T")
  if (is.matrix(mod$Phihat)) stop("Model was estimated using only one penalty parameter. Set cv=FALSE in sparseVAR.")

  n3 <- dim(mod$Phihat)[3]
  ics <- matrix(ncol = 3, nrow = n3)
  for (i in 1:n3){
    mod_tmp <- mod
    mod_tmp$Phihat <- mod$Phihat[, , i]
    mod_tmp$phi0hat <- mod$phi0hat[, , i]

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
    cat("\n\n#### Selected the following ####\n\n")
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
#' @param tt number of time periods
#' @param df degrees of freedom
.get_ic_vals <- function(res, tt, df){
  Omega <- cov(res)
  log_det_omega <- log(det(Omega))

  aic <- log_det_omega + 2/tt * df
  bic <- log_det_omega + log(tt)/tt * df
  hq <- log_det_omega + 2*log(log(tt))/tt * df

  out <- cbind(aic, bic, hq)

  out
}


