require(skewtDist) # rast
require(numDeriv) # jacobian, hessian
rm(list = ls())
source("ast-mle.R") # llast, llast_grad
source("helper-functions.R") # K, D, Dprime, L, R

set.seed(123)
n <- 1e4
mu <- 0.12
scale <- 0.6
alpha <- 0.6
nu1 <- 5
nu2 <- 6
pars <- c(mu, scale, alpha, nu1, nu2)
data <- rast(n, pars = pars)
arglist <- list(data = data,
                fixed_pars = rep(NA, 5))

llast(pars, arglist)

# jacobian
llast_grad(pars, arglist)
jacobian(llast, x = pars, arglist = arglist)

# hessian
hessian(llast, x = pars, arglist = arglist) / n

source("ast-infomat.R")
astInfoMat(pars = pars, method = "expected")
astInfoMat(pars = pars, data = data, method = "observed")




res <- astMLE(data = retSW[,1])


res$fitted_pars
astInfoMat(res$fitted_pars, data = retSW[,1])
astInfoMat(res$fitted_pars, data = retSW[,1], method = "observed")
