#' Plot the Cross Validation Error Curve for a Sparse VAR or VARX
#' @param fit Fitted VAR or VARX model.
#' returned by \code{\link{sparseVAR}} or \code{\link{sparseVARX}}.
#' @param ... Not currently used
#' @export
plot_cv <- function(fit, ...) {
  if ("bigtime.VAR" %in% class(fit)) {
    if (fit$selection != "cv") stop("No cross-validation was used in model estimation. Set selection='cv' in sparseVAR")
    # This copied mostly from ggb package
    m <- fit$MSFEcv
    se <- apply(fit$MSFE_all, 2, sd) / sqrt(nrow(fit$MSFE_all))
    yrang = range(c(m - se, m + se), na.rm = TRUE)
    par(xpd = FALSE)
    graphics::plot(fit$lambdas, m, xlab = "lambda",
                   ylab = "Cross-validation Error",
                   pch = 19, col = 2, ylim = yrang, log = "x", ...)
    error_bars(fit$lambdas, m - se, m + se, width = 0.01,
               col = "darkgrey")
    # graphics::abline(v = c(fit$lambda_opt, fit$lambda_SEopt), lty = 3)
    graphics::abline(v = c(fit$lambda_opt), lty = 3, lwd = 2)
    graphics::abline(v = c(fit$lambda_SEopt), lty = 3, lwd = 2, col = 2)
    graphics::legend("top", xpd = TRUE,
                     inset = c(0.05, -0.25),
                     legend = c("Opt. lambda", "SE Opt. lambda"),
                     col = c(1, 2), lty = c(3, 3), bty = "n", horiz = TRUE, lwd = c(2, 2))
  }
  else if ("bigtime.VARX" %in% class(fit)) {
    image_cv(fit)
  } else stop("Unsupported type.")
  invisible()
}

error_bars <- function (x, upper, lower, width = 0.02, ...) {
  # this function came from hierNet package
  xlim <- range(x)
  barw <- diff(xlim) * width
  graphics::segments(x, upper, x, lower, ...)
  graphics::segments(x - barw, upper, x + barw, upper, ...)
  graphics::segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

image_cv <- function(fit) {
  o <- order(unique(fit$lambdaB))
  graphics::image(x = unique(fit$lambdaB)[o],
        y = unique(fit$lambdaPhi),
        z = matrix(fit$MSFEcv, nrow = length(o))[o, ],
        xlab = "lambda_B", ylab="lambda_Phi", col = heat.colors(40))
  points(fit$lambdaB_SEopt, fit$lambdaPhi_SEopt, pch = 19)
}
