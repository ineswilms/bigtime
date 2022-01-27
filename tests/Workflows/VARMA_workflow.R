library(bigtime)
rm(list = ls())
data(varma.example)
Y <- Y.varma

########## Checking functionality via workflows ##################################################
mod <- sparseVARMA(Y) # This should return a tensor for phase II results and give a standardisation warning
stopifnot(length(dim(mod$Phihat))==3) # This should be a tensor
dim(mod$Phihat)
stopifnot(length(dim(mod$VARPhihat))==2) # This should be a matrix since first stage always does selection
dim(mod$VARPhihat)
stopifnot(mod$VARselection=="cv") # should be cv
stopifnot(mod$VARMAselection=="none") # should be none
stopifnot(class(mod)=="bigtime.VARMA") # Should be bigtime.VARMA

mod <- sparseVARMA(scale(Y)) # Same results as above but no standardisation warnings
dim(mod$Phihat) # Should be tensor
dim(mod$VARPhihat) # Should be matrix
res <- residuals(mod) # This should be a tensor
dim(res)
fit <- fitted(mod)
dim(fit) == dim(res)
dif <- sapply(1:dim(fit)[3], function(i) fit[,,i] + res[,,i] - mod$Y[-(1:max(mod$VARMAp, mod$VARMAq)), ], simplify = 'array')
max(abs(dif)) # This should be pretty much zero
try({lagmatrix(mod)}) # Should give error since no seleciton was used
try({plot_cv(mod)}) # Should give an error
try({diagnostics_plot(mod)}) # This should also give an error
dir_fcst <- directforecast(mod, h = 1)
dim(dir_fcst) # This should be a tensor

mod <- sparseVARMA(scale(Y), VARselection = "bic")
mod$VARselection
dim(mod$Phihat) # Should still be tensor
dim(mod$VARPhihat) # Should still be matrix

mod <- sparseVARMA(scale(Y), VARselection = "bic", VARMAselection = "cv")
mod$VARMAselection
dim(mod$Phihat) # Should be matrix
dim(mod$VARPhihat) # Should be matrix
res <- residuals(mod)
dim(res) # Should be a matrix
fit <- fitted(mod)
dim(fit) == dim(res)
dim(fit) # Should be a mtrix
max(abs(fit + res - mod$Y[-(1:max(mod$VARMAp, mod$VARMAq)), ])) # Close to zero
lagmatrix(mod, returnplot = TRUE)
plot_cv(mod)
diagnostics_plot(mod)


mod <- sparseVARMA(scale(Y), VARselection = "hq", VARMAselection = "bic") # Should print two IC tables
mod$VARMAselection
mod$VARselection
dim(mod$Phihat) # Should be matrix
dim(mod$VARPhihat) # Should be matrix
dir_fcst <- directforecast(mod, h = 1)
length(dir_fcst) # This should be a vector
dir_fcst


mod_lasso <- sparseVARMA(Y, VARpen = "L1") # STD warning
dim(mod_lasso$Phihat) # Should be tensor
dim(mod_lasso$VARPhihat) # Should be matrix


mod <- sparseVARMA(scale(Y), VARpen = "L1", VARMApen = "L1")
dim(mod$Phihat) # Should be tensor
dim(mod$VARPhihat) # Should be matrix

mod <- sparseVARMA(scale(Y), VARpen = "L1", VARMApen = "L1", VARMAselection = "cv")
dim(mod$Phihat) # Should be matrix
dim(mod_lasso$VARPhihat) # Should be matrix


mod <- sparseVARMA(scale(Y[, 1])) # Univariate case
dim(mod$Phihat) # Should be tensor
mod <- sparseVARMA(scale(Y[, 1]), VARMAselection = "cv")
dim(mod$Phihat) # Should be matrix
mod <- sparseVARMA(scale(Y[, 1]), VARMAselection = "bic")
dim(mod$Phihat) # Should be matrix
res <- residuals(mod)
dim(res) # Should be matrix
fit <- fitted(mod)
dim(fit) == dim(res)
max(abs(fit + res - mod$Y[-(1:max(mod$VARMAp, mod$VARMAq)), ])) # Close to zero
try({lagmatrix(mod, returnplot = TRUE)}) # Should give an error
try({plot_cv(mod)}) # Should give an error because not yet implemented
diagnostics_plot(mod)
dir_fcst <- directforecast(mod, h = 1)
dir_fcst

# Following test is from issue 4 on github
# The original version was not including the residual call
# but the error was in the residual funciton. Since the estimation
# procedure does not actually need the residual call, it was removed
# and instead the call was made explicit in the test.
testdata <- matrix(rnorm(300, 1:300), 100)
varmamod <- sparseVARMA(Y = scale(testdata), h=2)
res <- residuals(varmamod)


