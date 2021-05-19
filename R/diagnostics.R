
#' Gives the residuals for VAR models estimated using sparseVAR
#'
#' @param object model estimated using sparseVAR
#' @param ... Not currently used
#' @export
#' @return Returns a matrix of residuals.
#' @examples
#' dat <- simVAR(200, 2, 5, decay = 0.001, seed = 6150533)
#' mod <- sparseVAR(dat$Y)
#' f <- fitted(mod)
#' res <- resid(mod)
residuals.bigtime.VAR <- function(object, ...){
  mod <- object
  fit <- fitted.bigtime.VAR(mod, ...)
  s <- dim(mod$Y)[1] - dim(fit)[1]
  res <- mod$Y[-(1:s), ] - fit
  colnames(res) <- colnames(mod$Y)
  res
}

#' Gives the residuals for VAR models estimated using sparseVAR
#'
#' @param object model estimated using sparseVAR
#' @param ... Not currently used
#' @export
#' @return Returns a matrix of residuals.
#' @examples
#' dat <- bigtime::simVAR(200, 2, 5, decay = 0.001, seed = 6150533)
#' mod <- bigtime::sparseVAR(dat$Y)
#' f <- fitted(mod)
#' res <- resid(mod)
resid.bigtime.VAR <- function(object, ...){
  residuals.bigtime.VAR(object, ...)
}


#' Gives the fitted values of a model estimated using sparseVAR
#'
#' @param object model estimated using sparseVAR
#' @param ... Not currently used
#' @export
#' @return Returns a matrix of fitted values
#' @examples
#' dat <- bigtime::simVAR(200, 2, 5, decay = 0.001, seed = 6150533)
#' mod <- bigtime::sparseVAR(dat$Y)
#' f <- fitted(mod)
#' res <- resid(mod)
#'
fitted.bigtime.VAR <- function(object, ...){
  mod <- object
  VARdata <- HVARmodel(Y=mod$Y, p=mod$p, h=mod$h)
  fit <- t(apply(VARdata$fullZ, 2, function(x) mod$Phihat%*%x + mod$phi0hat))
  colnames(fit) <- colnames(mod$Y)
  fit
}


#' Creates a diagnostic Plot
#'
#' @param mod VAR model estimated using sparseVAR
#' @param variable Variable to show. Either numeric (which column) or character
#' (variable name)
#' @param dates Optional Date vector.
#' @export
#' @return Returns a ggplot2 plot
#' @examples
#' dat <- simVAR(200, 2, 5, decay = 0.1, seed = 6150533,
#'                        sparsity_pattern = "hvar")
#' mod <- sparseVAR(dat$Y, h = 1)
#' diagnostics_plot(mod, variable = 1)
diagnostics_plot <- function(mod, variable = 1, dates=NULL){
  if(!"bigtime.VAR" %in% class(mod)) stop("Only implemented for VAR models")
  fit <- fitted.bigtime.VAR(mod)
  res <- residuals.bigtime.VAR(mod)
  s <- dim(mod$Y)[1] - dim(fit)[1]
  Y <- mod$Y[-c(1:s), ]
  if (is.numeric(variable)) {
    if (variable > ncol(fit)) stop("Data does not have ", variable, " columns, but only ", ncol(fit))
    if (variable <= 0) stop("variable must be a name or an index >= 1")
  }else if(is.character(variable)) {
    if (is.null(colnames(fit))) stop("Data does not have column names.")
    else if (!(variable %in% colnames(fit))) stop("Data has no variable called ", variable)
    variable <- which(colnames(fit) == variable)
  }
  if (is.null(colnames(fit))){
    warning("Data has no column names. Using defaults Y1, ...")
    colnames(fit) <- paste0("Y", 1:ncol(fit))
    colnames(res) <- paste0("Y", 1:ncol(fit))
    colnames(Y) <-  paste0("Y", 1:ncol(fit))
  }
  variable <- colnames(fit)[variable]
  colnames(fit) <- paste0(colnames(fit), "_Fitted")
  colnames(Y) <- paste0(colnames(Y), "_Observed")
  colnames(res) <- paste0(colnames(res), "_Residual")
  df <- as.data.frame(cbind(fit, res, Y))

  d <- 1:nrow(fit) + s
  if(!is.null(dates)) {
    if (length(dates) != nrow(fit)) stop("dates must be of length ", nrow(fit), "not ", length(dates))
    d <- as.Date(dates)
  }

  cols <- c("Fitted"="red",
            "Observed"="black",
            "Residual"="black")
  Series <- Time <- Variable <- yintercept <- value <- p <- NULL # needed because of non-standard evaluation
  hline <- data.frame(list(p = "Residual", yintercept = 0))
  df %>%
    dplyr::mutate(Time = d) %>%
    tidyr::pivot_longer(-Time,
                        names_to = c("Variable", "Series"),
                        names_sep = "_") %>%
    dplyr::mutate(p = ifelse(Series %in% c("Fitted", "Observed"),
                             paste("Fitted", variable), "Residual"),
                  Series = factor(Series, levels = c("Observed",
                                                     "Fitted",
                                                     "Residual"))) %>%
    dplyr::filter(Variable == variable) %>%
    ggplot2::ggplot() +
    ggplot2::geom_hline(data = hline,
                        ggplot2::aes(yintercept = yintercept),
                        color = "red", linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(Time, value, color = Series)) +
    ggplot2::facet_grid(rows = dplyr::vars(p), scale = "free_y") +
    ggplot2::ylab("") +
    ggplot2::labs(color = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "top",
      legend.margin = ggplot2::margin(b=-10),
      # strip.background = ggplot2::element_blank(),
      # strip.text = ggplot2::element_blank()
    ) +
    ggplot2::scale_color_manual(values = cols, breaks = c("Observed", "Fitted"))
}


