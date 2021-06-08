
rm(list = ls())
library(bigtime)
data(Y)
data(X)

########## Simple Estimation ##################################################
mod_hlag <- sparseVARX(Y = scale(Y), X = scale(X), selection = "none") # estimation using hlag and no selection
dim(mod_hlag$Phihat)
mod_cv <- sparseVARX(Y = Y, X = X, selection = "cv") # estimation using cv
dim(mod_cv$Phihat)
mod_bic <- sparseVARX(Y = Y, X = X, selection = "bic") # estimation using bic
dim(mod_bic$Phihat)
mod_aic <- sparseVARX(Y = Y, X = X, selection = "aic") # estimation using aic
dim(mod_aic$Phihat)
mod_hq <- sparseVARX(Y = Y, X = X, selection = "hq") # estimation using hq
dim(mod_hq$Phihat)


mod_lasso <- sparseVARX(Y = Y, X = X, VARXpen = "L1") # estimation using lasso
class(mod_hlag)

########## diagnostics ########################################################

# Getting residuals and fitted values
# If it works for hlag then it also works for lasso
res_hlag <- residuals(mod_hlag)
fit_hlag <- fitted(mod_hlag)
all.equal(fit_hlag + res_hlag, Y[-(1:max(mod_hlag$p, mod_hlag$s)), ]) # This should be true otherwise something is wrong

# lagmatrix
lagmat <- lagmatrix(mod_cv, returnplot = TRUE) # This should give two plots, so go back one plot

# Checking cv procedure
plot_cv(mod_hlag)

# diagnostics plot
diagnostics_plot(mod_hlag)


########## forecasting ########################################################

# directforecasting --> Recursive forecasting is not possible for VARX
dir_fcst <- directforecast(mod_hlag, h = 1)
dir_fcst

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
plot_cv(mod) # should give an error
diagnostics_plot(mod) # this should throw an error

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
plot_cv(mod)
diagnostics_plot(mod)
