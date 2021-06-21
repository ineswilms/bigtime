#' It is difficult to come up with unit tests for the statistical funcitons
#' included in bigtime. For these reasons, this document represents an
#' informal but yet unified way to test the package after larger updates.

library(bigtime)
rm(list = ls())
######### settings ############################################################
seed <- 6150533

######### simulating data (no sparsity) #######################################

k <- 2
p <- 2
periods <- 250
max_abs_eigval <- 0.8
e_dist <- function(n) rnorm(n, mean = 0, sd = 0.001)
sim_dat <- simVAR(periods = periods, k = k, p = p, e_dist = e_dist,
                  seed = seed, max_abs_eigval = max_abs_eigval)

#' Tests for the coefficient matrix exist. We hence only need to focus on
#' the statistical parts of the data simulation, like e_dist and
#' whether we truly have the VAR we were expecting

# Checking e_dist
e <- sim_dat$e_dist[1:k, ]
e_mu <- apply(e, 1, mean) # should be close to true mean
e_sd <- apply(e, 1, sd) # should be close to true sd
cat(paste("True/Est\t Mean:",  0, " / ", e_mu, "\t SD:", 0.001, " / ", e_sd, "\n"), sep = "")

# Checking the VAR visually
plot.ts(sim_dat$Y)

# Checking whether we can recover the VAR
require(vars)
mod <- vars::VAR(sim_dat$Y, p = p, type = "const")
true.coefs <- capture.output(print(cbind(sim_dat$coef_mat[1:k, ], sim_dat$const[1:k])))
est.coefs <- capture.output(print(vars::Bcoef(mod)))
cat("True", true.coefs, "Estimated", est.coefs, sep = "\n")

######### simulating data (lasso) #############################################

k <- 3
p <- 3
periods <- 250
max_abs_eigval <- 0.8
e_dist <- function(n) rnorm(n, mean = 0, sd = 0.001)
sim_dat <- simVAR(periods = periods, k = k, p = p, e_dist = e_dist,
                             seed = seed, max_abs_eigval = max_abs_eigval,
                             sparsity_pattern = "lasso",
                             sparsity_options = list(num_zeroes = 15))
# We can simply check the coefficient matrix on whether it has the right struct.
print(sim_dat$coef_mat[1:k, ])
print(sum(sim_dat$coef_mat[1:k, ] == 0))


######### simulating data (hvar) ##############################################

k <- 3
p <- 3
periods <- 250
max_abs_eigval <- 0.8
e_dist <- function(n) rnorm(n, mean = 0, sd = 0.001)
sim_dat <- simVAR(periods = periods, k = k, p = p, e_dist = e_dist,
                             seed = seed, max_abs_eigval = max_abs_eigval,
                             sparsity_pattern = "hvar",
                             sparsity_options = list(
                               zero_min = 0,
                               zero_max = 3,
                               zeroes_in_self = TRUE
                             ))


# Inspecting the coef_mat
print(sim_dat$coef_mat[1:k, ])


######### simulating univariate data ##########################################

k <- 1
p <- 3
periods <- 250
max_abs_eigval <- 0.8
e_dist <- function(n) rnorm(n, mean = 0, sd = 0.001)
sim_dat <- simVAR(periods = periods, k = k, p = p, e_dist = e_dist,
                  seed = seed, max_abs_eigval = max_abs_eigval,
                  sparsity_pattern = "hvar",
                  sparsity_options = list(
                    zero_min = 0,
                    zero_max = 3,
                    zeroes_in_self = TRUE
                  ))
sim_dat$coef_mat
plot.ts(sim_dat$Y)

######### Estimation (hvar) ###################################################

k <- 3
p <- 3
periods <- 250
max_abs_eigval <- 0.8
e_dist <- function(n) rnorm(n, mean = 0, sd = 0.001)
sim_dat <- simVAR(periods = periods, k = k, p = p, e_dist = e_dist,
                  seed = seed, max_abs_eigval = max_abs_eigval,
                  sparsity_pattern = "hvar",
                  sparsity_options = list(
                    zero_min = 0,
                    zero_max = 3,
                    zeroes_in_self = TRUE
                  ))

# Using no selection at all
hvar <- sparseVAR(scale(sim_dat$Y))
# hvar should not use any selection method.
# Thus selection should be none and phis should be a 3d array
hvar$selection
dim(hvar$Phihat)
dim(hvar$phi0hat)

hvar <- sparseVAR(sim_dat$Y) # This should give a standardisation warning

# Using cv selection
hvar <- sparseVAR(scale(sim_dat$Y), selection = "cv")
# selection field should be cv, lambda opt and lambda SE opt should be given
# and phis should be matrix/vector
hvar$selection
dim(hvar$Phihat)
length(hvar$phi0hat)
hvar$lambda_opt
hvar$lambda_SEopt

# using bic
hvar <- sparseVAR(scale(sim_dat$Y), selection = "bic")
# selection should be bic, and phis should be matrix/vector
# lambda opt should contain bic opt choice
# verbose is always true, so output should be shown comparing different lambdas
# and ic
hvar$selection
hvar$lambda_opt
hvar$lambda_SEopt # This should be NA
dim(hvar$Phihat)
length(hvar$phi0hat)

# using aic
hvar <- sparseVAR(scale(sim_dat$Y), selection = "aic")
# selection should be aic, and phis should be matrix/vector
# lambda opt should contain aic opt choice
# verbose is always true, so output should be shown comparing different lambdas
# and ic
hvar$selection
hvar$lambda_opt
hvar$lambda_SEopt # This should be NA
dim(hvar$Phihat)
length(hvar$phi0hat)

# using hq
hvar <- sparseVAR(scale(sim_dat$Y), selection = "hq")
# selection should be hq, and phis should be matrix/vector
# lambda opt should contain hq opt choice
# verbose is always true, so output should be shown comparing different lambdas
# and ic
hvar$selection
hvar$lambda_opt
hvar$lambda_SEopt # This should be NA
dim(hvar$Phihat)
length(hvar$phi0hat)

######### Estimation (lasso) ##################################################
# If above looks all good, then we can check whether all also works for lasso

# no selection method
lasso <- sparseVAR(scale(sim_dat$Y), selection = "none", VARpen = "L1")
lasso$selection
dim(lasso$Phihat)
dim(lasso$phi0hat)

# cv selection
lasso <- sparseVAR(scale(sim_dat$Y), selection = "cv", VARpen = "L1")
lasso$selection
lasso$lambda_opt
lasso$lambda_SEopt
dim(lasso$Phihat)
length(lasso$phi0hat)


# only checking bic selection
lasso <- sparseVAR(scale(sim_dat$Y), selection = "bic", VARpen = "L1")
lasso$selection
lasso$lambda_opt
lasso$lambda_SEopt
dim(lasso$Phihat)
length(lasso$phi0hat)

######### VAR diagnostics #####################################################
k <- 3
p <- 3
periods <- 250
max_abs_eigval <- 0.8
e_dist <- function(n) rnorm(n, mean = 0, sd = 0.001)
sim_dat <- simVAR(periods = periods, k = k, p = p, e_dist = e_dist,
                  seed = seed, max_abs_eigval = max_abs_eigval,
                  sparsity_pattern = "hvar",
                  sparsity_options = list(
                    zero_min = 0,
                    zero_max = 3,
                    zeroes_in_self = TRUE
                  ))


hvar <- sparseVAR(scale(sim_dat$Y), selection = "none")
try({plot_cv(hvar)}) # Should give an error because model did not use cv
try({diagnostics_plot(hvar)}) # This should not work and should give an error
res <- residuals(hvar)
dim(res) # Should be a 3d array
fit <- fitted(hvar)
dim(fit) # Should have same dimensions as res
try({lagmat <- lagmatrix(hvar, returnplot = TRUE)}) # This should not work

hvar <- sparseVAR(scale(sim_dat$Y), selection = "cv")
plot_cv(hvar) # This should work
diagnostics_plot(hvar) # This should work and plot for the first series
res <- resid(hvar) # Can also be used instead of residuals
dim(res) # Should be a matrix. If it works for this it will also work for all IC selections
fit <- fitted(hvar)
dim(fit) # Should have same dimensions as res
lagmat <- lagmatrix(hvar, returnplot = TRUE) # This should work just fine.

hvar <- sparseVAR(scale(sim_dat$Y), selection = "bic")
try({plot_cv(hvar)}) # Should not work
diagnostics_plot(hvar, variable = 2) # This should work and plot for the second series
lagmat <- lagmatrix(hvar, returnplot = TRUE) # This should work just fine. If it works here, it also works for the followign

hvar <- sparseVAR(scale(sim_dat$Y), selection = "aic")
try({plot_cv(hvar)}) # Should not work
diagnostics_plot(hvar, variable = 3) # This should work and plot for the third series

hvar <- sparseVAR(scale(sim_dat$Y), selection = "hq")
try({plot_cv(hvar)}) # Should not work
diagnostics_plot(hvar, variable = 'Y2') # This should work and plot for the second series
try({diagnostics_plot(hvar, variable = 4)}) # This should not work since only 3 series are estimated
try({diagnostics_plot(hvar, variable = "Y4")}) # This should not work since only 3 series are estimated

# It should also all work if we use lasso instead
hvar <- sparseVAR(scale(sim_dat$Y), selection = "cv", VARpen = "L1")
plot_cv(hvar) # This should work
diagnostics_plot(hvar) # This should work and plot for the first series
res <- resid(hvar) # Can also be used instead of residuals
dim(res) # Should be a matrix. If it works for this it will also work for all IC selections
fit <- fitted(hvar)
dim(fit) # Should have same dimensions as res
lagmat <- lagmatrix(hvar, returnplot = TRUE) # This should work just fine.

######### VAR forecasting #####################################################
k <- 3
p <- 3
periods <- 250
max_abs_eigval <- 0.8
e_dist <- function(n) rnorm(n, mean = 0, sd = 0.001)
sim_dat <- simVAR(periods = periods, k = k, p = p, e_dist = e_dist,
                  seed = seed, max_abs_eigval = max_abs_eigval,
                  sparsity_pattern = "hvar",
                  sparsity_options = list(
                    zero_min = 0,
                    zero_max = 3,
                    zeroes_in_self = TRUE
                  ))

# doing it for hvar
mod <- sparseVAR(scale(sim_dat$Y), selection = "none")
dir_fcst <- directforecast(mod)
dim(dir_fcst) # This should be a 3d array
rec_fcst <- recursiveforecast(mod, h = 5) # message should be given saying all lambdas are used for forecasting
dim(rec_fcst$fcst) # This should be a 3d array
plot(rec_fcst) # This should plot a ribbon covering the min to max forecast at each horizon

# Same should work if lasso is done
mod <- sparseVAR(scale(sim_dat$Y), selection = "none", VARpen = "L1")
dir_fcst <- directforecast(mod)
dim(dir_fcst) # This should be a 3d array
rec_fcst <- recursiveforecast(mod, h = 5) # message should be given saying all lambdas are used for forecasting
dim(rec_fcst$fcst) # This should be a 3d array
plot(rec_fcst) # This should plot a ribbon covering the min to max forecast at each horizon

# selecting using cv and any ic should behave same and thus only cv is tested
mod <- sparseVAR(scale(sim_dat$Y), selection = "cv")
dir_fcst <- directforecast(mod)
length(dir_fcst) # This should be a vector
rec_fcst <- recursiveforecast(mod, h = 5)
dim(rec_fcst$fcst) # This should be a matrix
plot(rec_fcst) # Should plot a line since only one lambda is used


###### Testing univariate functionality #######################################
k <- 1
p <- 3
periods <- 250
max_abs_eigval <- 0.8
e_dist <- function(n) rnorm(n, mean = 0, sd = 0.001)
sim_dat <- simVAR(periods = periods, k = k, p = p, e_dist = e_dist,
                  seed = seed, max_abs_eigval = max_abs_eigval,
                  sparsity_pattern = "hvar",
                  sparsity_options = list(
                    zero_min = 0,
                    zero_max = 3,
                    zeroes_in_self = TRUE))


## Estimation using none
mod <- sparseVAR(scale(sim_dat$Y), selection = "none")
dim(mod$phi0hat)
dim(mod$Phihat)
fit <- fitted(mod)
dim(fit)
res <- residuals(mod)
dim(res)
max(abs(fit + res - mod$Y[-(1:mod$p), ])) # this should be zero
try({plot_cv(mod)}) # this should throw an error
try({diagnostics_plot(mod)}) # this should also throw an error

## Estimation using cv
mod <- sparseVAR(scale(sim_dat$Y), selection = "cv")
dim(mod$phi0hat)
length(mod$phi0hat)
dim(mod$Phihat)
length(mod$Phihat)
fit <- fitted(mod)
res <- residuals(mod)
max(abs(fit + res - mod$Y[-(1:mod$p), ]))
plot_cv(mod)
diagnostics_plot(mod, variable = 1)


## Estmimation using bic, aic, hq
mod <- sparseVAR(scale(sim_dat$Y), selection = "bic")
mod <- sparseVAR(scale(sim_dat$Y), selection = "aic")
mod <- sparseVAR(scale(sim_dat$Y), selection = "hq")
mod$selection
mod$lambda_opt
try({plot_cv(mod)}) # This should throw an error
diagnostics_plot(mod)
