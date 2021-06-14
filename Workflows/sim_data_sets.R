library(bigtime)
rm(list = ls())


# VAR data simulation
?bigtime::simVAR
periods <- 200
k <- 5
p <- 8
seed <- 6150533
sparsity_pattern <- "hvar"
sparsity_options = list(
  zero_min = 3,
  zero_max = 8,
  zeroes_in_self = TRUE
)

VAR <- bigtime::simVAR(periods, k, p,
                       sparsity_pattern = sparsity_pattern,
                       sparsity_options = sparsity_options,
                       seed = seed)
summary(VAR)
mod <- bigtime::sparseVAR(scale(VAR$Y), selection = "bic")
lagmatrix(mod, returnplot = TRUE)
Y <- VAR$Y
save(Y, file = "./data/var.example.rda")


# VARX data simulation
# 1. Simulate Exogenous Data
# 2. create edist as composite of error and exogenous data
# 3. simulate VAR with custom e_dist
set.seed(6150533)
periods <- 200
burnin <- 200
k.x <- 3
p.x <- 1
k <- 3
p <- 4
seed <- 6150533
sparsity_pattern <- "hvar"
sparsity_options = list(
  zero_min = 0,
  zero_max = 1,
  zeroes_in_self = TRUE
)
X.sim <- simVAR(periods+burnin+1, k.x, p.x,
                sparsity_pattern = sparsity_pattern,
                sparsity_options = sparsity_options,
                seed = seed)
summary(X.sim)
X <- lag(X.sim$Y)[-1, ]
coef_mat <- create_rand_coef_mat(k.x, p.x,
                                 sparsity_pattern = sparsity_pattern,
                                 sparsity_options = sparsity_options,
                                 max_abs_eigval = 0.95)
e <- rnorm((periods+burnin)*k, sd = 0.1)
e <- matrix(e, nrow = k)
edist <- coef_mat%*%t(X) + e
sparsity_options = list(
  zero_min = 0,
  zero_max = 4,
  zeroes_in_self = TRUE
)
Y.sim <- simVAR(periods, k, p,
                burnin = burnin,
                sparsity_pattern = sparsity_pattern,
                sparsity_options = sparsity_options,
                seed = seed)
summary(Y.sim)
Y <- Y.sim$Y
X <- X.sim$Y[-(1:(burnin+1)), ]
mod <- sparseVARX(Y = scale(Y), X = scale(X), selection = "cv")
plot_cv(mod)
lagmatrix(mod, returnplot = TRUE)

colnames(X) <- paste0("X", 1:ncol(X))
save(X, Y, file = "./data/varx.example.rda")


# VARMA data simulation
# 1. Simulate Errors
# 2. embed errors in right dimension and form e_dist
# 3. simulate VAR with custom e_dist
set.seed(14315)
k <- 3
p <- 4
p.e <- 1
burnin <- 200
periods <- 200
u.lag.coefs <- c(0.5, 0.8, 0.1)
sparsity_options = list(
  zero_min = 0,
  zero_max = 4,
  zeroes_in_self = TRUE
)

e <- rnorm((periods+burnin+1)*k, sd = 0.1)
e <- matrix(e, nrow = k)
u <- lapply(1:k, function(i) embed(e[i, ], dimension = p.e + 1))
e_dist <- lapply(1:k, function(i) rowSums(t(t(u[[i]]) * c(1, u.lag.coefs[i]))))
e_dist <- do.call(rbind, e_dist)

Y.sim <- simVAR(periods, k, p, e_dist = e_dist,
                sparsity_pattern = "hvar",
                sparsity_options = sparsity_options,
                burnin = burnin,
                seed = 6150533)
summary(Y.sim)
Y <- Y.sim$Y
save(Y, file = "./data/varma.example.rda")


mod <- sparseVARMA(scale(Y), VARMAselection = "cv")
lagmatrix(mod, returnplot = TRUE)
diagnostics_plot(mod)
plot_cv(mod)






