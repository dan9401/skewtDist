require(skewtDist) # rast
require(numDeriv) # jacobian, hessian
rm(list = ls())
source("ast-mle.R") # llast, llast_grad
source("helper-functions.R") # K, D, Dprime, L, R


set.seed(22)
n <- 1e7
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

# Dprime is exclusively used for observed method
# munu1, munu2
# sigma2, sigmanu1, sigmanu2
# alphanu1, alphanu2
# nu12 nu22
astInfoMat(pars = pars, data = data, method = "observed")

y <- data
y1 <- y[y <= mu]
y2 <- y[y > mu]

S_alpha2 <- length(y1)/length(y) * ((nu1 + 1) / alpha)^2 * mean( ( 1 - 1/L(pars, y1) )^2 ) +
  length(y2)/length(y) * ((nu2 + 1) / (1 - alpha))^2 * mean((1 - 1/R(pars, y2))^2)

S_alpha21 <- (nu1+1)/alpha^2 * mean(1 + 1/L(pars, y1) - 2/L(pars, y1)^2) +
  (nu2+1)/(1-alpha)^2 * mean(1 + 1/R(pars, y2) - 2/R(pars, y2)^2)


res <- astMLE(data = retSW[,1])


res$fitted_pars
astInfoMat(res$fitted_pars, data = retSW[,1])
astInfoMat(res$fitted_pars, data = retSW[,1], method = "observed")
