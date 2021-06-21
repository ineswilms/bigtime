rm(list = ls())
library(bigtime)
data(varx.example)
Y <- Y.varx
X <- X.varx

########## Simple Estimation ##################################################
mod_hlag <- sparseVARX(Y = scale(Y), X = scale(X), selection = "none") # estimation using hlag and no selection
dim(mod_hlag$Phihat) # should be tensor
try({plot_cv(mod_hlag)}) # Should throw error
try({lagmatrix(mod_hlag, returnplot = TRUE)}) # Should throw error
try({diagnostics_plot(mod_hlag)}) # Should throw error
res <- residuals(mod_hlag)
dim(res) # Should be tensor
fit <- fitted(mod_hlag)
dim(fit) # Should be tensor
max(abs(sapply(1:dim(fit)[3], function(i) fit[,,i] + res[,,i] - mod_hlag$Y[-c(1:max(mod_hlag$p, mod_hlag$s)), ], simplify = 'array'))) # Should be close to zero
dir_fcst <- directforecast(mod_hlag, h = 1)
dim(dir_fcst) # Should be a tensor

# Should give standardisation warnings
mod_cv <- sparseVARX(Y = Y, X = X, selection = "cv") # estimation using cv
dim(mod_cv$Phihat) # Should be matrix
plot_cv(mod_cv)
lagmatrix(mod_cv, returnplot = TRUE)
diagnostics_plot(mod_cv)
res <- residuals(mod_cv)
dim(res) # Should be matrix
fit <- fitted(mod_cv)
dim(fit) # Should be matrix
max(abs(fit+res-mod_cv$Y[-c(1:max(mod_cv$p, mod_cv$s)), ])) # Should be close to zero
dir_fcst <- directforecast(mod_cv, h = 1)
length(dir_fcst) # Should be vector


mod_bic <- sparseVARX(Y = Y, X = X, selection = "bic") # estimation using bic
dim(mod_bic$Phihat)
try({plot_cv(mod_bic)}) # should throw error
mod_aic <- sparseVARX(Y = Y, X = X, selection = "aic") # estimation using aic
dim(mod_aic$Phihat)
mod_hq <- sparseVARX(Y = Y, X = X, selection = "hq") # estimation using hq
dim(mod_hq$Phihat)

mod_lasso <- sparseVARX(Y = Y, X = X, VARXpen = "L1") # estimation using lasso
class(mod_hlag)

########## univariate case ###################################################
Y <- Y[, 1, drop = FALSE]
X <- X[, 1, drop = FALSE]

mod <- sparseVARX(Y = scale(Y), X = scale(X), selection = "none")
dim(mod$phi0hat)
fit <- fitted(mod)
dim(fit)
res <- residuals(mod)
dim(res)
max(abs(fit + res - mod$Y[-(1:max(mod$p, mod$s)), ])) # should be  close to zero
try({plot_cv(mod)}) # should give an error
try({diagnostics_plot(mod)}) # this should throw an error

# CV estimation
mod <- sparseVARX(Y = scale(Y), X = scale(X), selection = "cv")
dim(mod$Phihat)
fit <- fitted(mod)
dim(fit)
res <- residuals(mod)
dim(res)
max(abs(fit + res - mod$Y[-(1:max(mod$p, mod$s)), ]))
plot_cv(mod)
diagnostics_plot(mod)

# estimation using IC
mod <- sparseVARX(Y = scale(Y), X = scale(X), selection = "bic")
mod <- sparseVARX(Y = scale(Y), X = scale(X), selection = "aic")
mod <- sparseVARX(Y = scale(Y), X = scale(X), selection = "hq")

fit <- fitted(mod)
res <- residuals(mod)
max(abs(fit + res - mod$Y[-c(1:max(mod$p, mod$s)), ]))
try({plot_cv(mod)}) # Should throw an error because no cv was done
diagnostics_plot(mod)
