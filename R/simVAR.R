#' Simulates a VAR(p) with various sparsity patterns
#'
#' @param periods Scalar indicating the desired time series length
#' @param k Number of time series
#' @param p Maximum lag number. In case of \code{sparsity_patter="none"} this will be
#' the actual number of lags for all variables
#' @param coef_mat Coefficient matrix in companion form. If not provided,
#' one will be simulated
#' @param const Constant term of VAR. Default is zero. Must be either a scalar,
#' in which case it will be broadcasted to a k-vector, or a k-vector
#' @param e_dist Either a function taking argument n indicating the number of
#' variables in the system, or a matrix of dimensions k x (periods+burnin)
#' @param init_y Initial values. Defaults to zero. Expects either a scalar or
#' a vector of length (k*p)
#' @param max_abs_eigval Maximum allowed eigenvalue of companion matrix.
#' Only applicable if coefficient matrix is being simulated
#' @param burnin Number of time points to be used for burnin
#' @param sparsity_pattern The sparsity pattern that should be simulated.
#' Options are: \code{"none"} for a dense VAR, \code{"lasso"} (or \code{"L1"})
#' for a VAR with random zeroes,
#' and \code{"hvar"} (or \code{"HLag"}) for an elementwise hierarchical sparsity pattern
#' @param sparsity_options Named list of additional options for
#' when sparsity pattern is lasso (L1) or hvar (HLag). For lasso (L1) the option \code{num_zero}
#' determines the number of zeros. For hvar (HLag), the options \code{zero_min} (\code{zero_max})
#' give the minimum (maximum) of zeroes for each variable in each equation,
#' and the option \code{zeroes_in_self} (boolean) determines if any of the
#' coefficients of a variable on itself should be zero.
#' @param decay How much smaller should parameters for later lags be. The
#' smaller, the larger will early parameters be w.r.t. later ones.
#' @param seed Seed to be used for the simulation
#' @param ... Additional arguments passed to \code{e_dist}
#' @export
#' @return Returns an object of S3 class \code{bigtime.simVAR} containing the following
#' \item{Y}{Simulated Data}
#' \item{periods}{Time series length}
#' \item{k}{Number of endogenous variables}
#' \item{p}{Maximum lag length; effective lag length might be shorter due to sparsity patterns}
#' \item{coef_mat}{Companion form of the coefficient matrix. Will be of
#' dimensions (\code{k}\code{p})x(\code{k}\code{p}). First \code{k} rows correspond
#' to the actual coefficient matrix.}
#' \item{is_coef_mat_simulated}{\code{TRUE} if the \code{coef_mat} was simulated, \code{FALSE} if
#' it was user provided}
#' \item{const}{Constant term}
#' \item{e_dist}{Errors used in the construction of the data}
#' \item{init_y}{Initial conditions}
#' \item{max_abs_eigval}{Maximum eigenvalue to which the companion matrix
#' was constraint}
#' \item{burnin}{Burnin period used}
#' \item{sparsity_pattern}{Sparsity pattern used}
#' \item{sparsity_options}{Extra options for the sparsity patterns used}
#' \item{seed}{Seed used for the simulation}
simVAR <- function(periods, k, p, coef_mat = NULL, const = rep(0, k), e_dist = rnorm,
                   init_y = rep(0, k*p), max_abs_eigval = 0.8, burnin = periods,
                   sparsity_pattern = c("none", "lasso", "L1", "hvar", "HLag"),
                   sparsity_options = NULL, decay = 1/p,
                   seed = NULL,
                   ...){

  ##################################
  #### Check and Prep ##############
  ##################################

  if (is.null(seed)) warning("No seed is being used. Replication might be infeasible")
  if (!is.null(seed)) set.seed(seed)

  # Checking all variables
  if (periods <= 0) stop("Must simulate at least one period: periods > 0")
  if (k < 1) stop("Must have at least one variable: k >= 1")
  if (!is.function(e_dist) & !is.matrix(e_dist)) stop("e_dist must either be a function or a matrix of dimensions k x (periods + burnin)")
  if (is.matrix(e_dist)){
    if ((ncol(e_dist) < periods | nrow(e_dist) < k)) stop("e_dist must be a function or a matrix of dimensions k x (periods + burnin)")
  }
  if (!is.null(coef_mat)) {
    if (k*p != ncol(coef_mat)) stop("coef_mat must have k*p columns not ", ncol(coef_mat))
  }
  if (!(length(const) == 1 | length(const) == k)) stop("const must either be of length one or of length k, not of length", length(const))
  if (p < 1) stop("Must include at least one lag: p >= 1")
  if (!(length(init_y) == 1 | length(init_y) == k*p)) stop("init_y must be of length one or of length k*p, not of length", length(init_y))
  if (max_abs_eigval < 0 | max_abs_eigval >=1) stop("Maximum eigenvalue must be strictly between zero and one, not ", max_abs_eigval)
  if (burnin < 0) stop("burnin must be non-negative: burnin >= 0")
  if (burnin == 0) warning("Not burning in the series is not recommended")
  if (!is.null(coef_mat)) warning("Will not use provided sparsity_pattern since coef_mat was given")
  if (decay > 1 | decay <=0) stop("Decay parameter must be between 0 and 1 (including)")

  # Which sparsity pattern was chosen?
  sparsity_pattern <- match.arg(sparsity_pattern)
  if (sparsity_pattern == "L1") sparsity_pattern <- "lasso"
  if (sparsity_pattern == "HLag") sparsity_pattern <- "hvar"

  # Creating coef_mat if needed
  is_coef_mat_simulated <- is.null(coef_mat)
  if (is.null(coef_mat)) coef_mat <- create_rand_coef_mat(k, p, max_abs_eigval, sparsity_pattern, sparsity_options, decay = decay)

  # Getting error terms
  if (is.function(e_dist)) e_dist <- do.call(cbind, lapply(1:(periods+burnin), function(x) e_dist(n = k, ...)))
  e_dist <- rbind(e_dist, matrix(0, nrow = k*(p-1), ncol = ncol(e_dist)))

  # constructing the constent
  if (length(const) == 1) const <- rep(const, k)
  const <- c(const, rep(0, k*(p-1)))

  # constructing the intial values
  if (length(init_y) == 1) init_y <- rep(init_y, k*p)

  ######################################
  #### END Check and Prep ##############
  ######################################


  Y <- simVAR_cpp(periods, k, p, coef_mat, const, e_dist, init_y, burnin)
  colnames(Y) <- paste0("Y", 1:k)

  sim_data <- list(
    Y = Y,
    periods = periods,
    k = k,
    p = p,
    coef_mat = coef_mat,
    is_coef_mat_simulated = is_coef_mat_simulated,
    const = const,
    e_dist = e_dist,
    init_y = init_y,
    max_abs_eigval = max_abs_eigval,
    burnin = burnin,
    sparsity_pattern = sparsity_pattern,
    sparsity_options = sparsity_options,
    seed = seed
  )
  class(sim_data) <- "bigtime.simVAR"
  sim_data
}


#' Creates a random coefficient matrix
#'
#' @param k Number of time series
#' @param p Number of lags
#' @param max_abs_eigval if < 1, then the VAR will be stable
#' @param sparsity_pattern The sparsity pattern that should be simulated.
#' Options are: \code{"none"} for a dense VAR, \code{"lasso"} for a VAR with random zeroes,
#' and \code{"hvar"} for an elementwise hierarchical sparsity pattern
#' @param sparsity_options Named list of additional options for
#' when sparsity pattern is lasso or hvar. For lasso the option \code{num_zero}
#' determines the number of zeros. For hvar, the options \code{zero_min} (\code{zero_max})
#' give the minimum (maximum) of zeroes for each variable in each equation,
#' and the option \code{zeroes_in_self} (boolean) determines if any of the
#' coefficients of a variable on itself should be zero.
#' @param decay How fast should coefficients shrink when the lag increases.
#' @param ... Not currently used
#' @export
#' @return Returns a coefficient matrix in companion form of dimension \code{kp}x\code{kp}.
create_rand_coef_mat <- function(k, p,
                                 max_abs_eigval = 0.8,
                                 sparsity_pattern = c("none", "lasso", "hvar"),
                                 sparsity_options = NULL,
                                 decay = 0.5,
                                 ...){

  ##################################
  #### Check and Prep ##############
  ##################################


  if (max_abs_eigval <= 0) stop("max_abs_eigval must be strictly postive, max_abs_eigval > 0")
  if (k < 1) stop("Model must have at least one variable: k>=1")
  if (p < 1) stop("Model must contain at least one lag: p > 1")

  sparsity_pattern <- match.arg(sparsity_pattern)
  # phi_vals <- runif(k*p*k, min = -100, max = 100)
  phi_vals <- rnorm(k*p*k, sd = 10)
  phis <- matrix(phi_vals, nrow = k, ncol = p*k)
  scale <- decay^(seq(0, p-1, by = 1))
  scale <- rep(scale, each = k)
  phis <- t(apply(phis, 1, function(x) x*scale))

  if (sparsity_pattern == "lasso"){
    # Some phi_vals will randomly be set to zero
    num_zero <- sparsity_options$num_zero
    if(is.null(num_zero)) {
      warning("Number of coefficients to be zero (num_zero) was not specified in sparsity_options. Defaulting to ", k*p)
      num_zero <- k*p
    }
    set_zero <- sample(1:(k*k*p), num_zero, replace = FALSE)
    phis[set_zero] <- 0
  }
  else if (sparsity_pattern == "hvar"){
    zero_min <- sparsity_options$zero_min
    zero_max <- sparsity_options$zero_max
    zeroes_in_self <- sparsity_options$zeroes_in_self
    if (is.null(zero_min)){
      warning("Minimum of zero was not given in sparsity_options. Defaulting to 0")
      zero_min <- 0
    }
    if (is.null(zero_max)){
      warning("Maximum of zero was not given in sparsity_options. Defaulting to ", floor(p/2))
      zero_max <- floor(p/2)
    }
    if (is.null(zeroes_in_self)){
      warning("Not specified whether there should be zero in self lags. Defaulting to true")
      zeroes_in_self <- TRUE
    }
    if(zero_max < zero_min) stop("zero_min must be smaller-equal than zero_max")
    if (zero_min > p | zero_max > p) stop("zero_min and zero_max must be at most p")

    for (k1 in 1:k){
      for (k2 in 1:k){
        if (k2 == k1 & !zeroes_in_self) next
        zeros <- sample(zero_min:zero_max, 1)
        if (zeros == 0) next
        set_zero <- seq(k2, k*p, k)[(p+1-zeros):p]
        phis[k1, set_zero] <- 0
      }
    }
  }

  ######################################
  #### END Check and Prep ##############
  ######################################


  I <- diag(1, nrow = k*(p-1), k*(p-1))
  O <- matrix(0, nrow = k*(p-1), ncol = k)
  coef_mat <- rbind(phis, cbind(I, O))
  while (max(abs(eigen(coef_mat)$values))>max_abs_eigval) {
    phis <- phis*0.99
    coef_mat <- rbind(phis, cbind(I, O))
  }
  rownames(coef_mat) <- c(paste0("Y", 1:k), rep("", dim(coef_mat)[[1]]-k))
  colnames(coef_mat) <- paste0(rownames(coef_mat)[1:k], ".L", rep(1:p, each = k))

  coef_mat
}

#' Plots a simulated VAR
#' @param x Simulated data of class \code{bigtime.simVAR} obtained
#' from the \code{simVAR}} function
#' @param ... Not currently used
#' @export
#' @return Returns a ggplot2 plot
plot.bigtime.simVAR <- function(x, ...){
  Time <- SimData <- Series <- NULL # Needed because of non-standard evaluation

  sim_data <- x
  Y <- sim_data$Y
  # if (!require(tidyverse, quietly = TRUE)) stop("tidyverse must be installed")
  dplyr::as_tibble(Y) %>%
    dplyr::mutate(Time = 1:dplyr::n()) %>%
    tidyr::pivot_longer(-Time, names_to = "Series", values_to = "SimData") %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(Time, SimData)) +
    ggplot2::facet_grid(rows = dplyr::vars(Series)) +
    ggplot2::ylab("") +
    ggplot2::theme_bw()
}


#' Gives a small summary of a VAR simulation
#' @param object Simulated data of class \code{bigtime.simVAR} obtained
#' from the \code\link{simVAR}} function
#' @param plot Should the VAR be plotted. Default is \code{TRUE}
#' @param ... Not currently used
#' @export
#' @return If `plot=TRUE`, then a ggplot2 plot will be returned
summary.bigtime.simVAR <- function(object, plot = TRUE, ...){
  sim_data <- object
  cat("#### General Information #### \n\n")
  cat("Seed                                       ", sim_data$seed, "\n")
  cat("Time series length                         ", sim_data$periods, "\n")
  cat("Burnin                                     ", sim_data$burnin, "\n")
  cat("Variables Simulated                        ", sim_data$k, "\n")
  cat("Number of Lags                             ", sim_data$p, "\n")
  cat("Coefficients were randomly created?        ", sim_data$is_coef_mat_simulated, "\n")
  cat("Maximum Eigenvalue of Companion Matrix     ", sim_data$max_abs_eigval, "\n")
  cat("Sparsity Pattern                           ", sim_data$sparsity_pattern, "\n")

  cat("\n\n#### Sparsity Options #### \n\n")
  print(sim_data$sparsity_options)

  cat("\n\n#### Coefficient Matrix #### \n\n")
  print(t(sim_data$coef_mat[1:sim_data$k, ]))

  if(plot) plot(sim_data)
}






