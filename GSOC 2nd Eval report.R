rm(list = ls())
devtools::install_github("dan9401/st")
library(st)
library(xts)

#
pars <- c(0.12, 0.6, 0.6, 3, 5)
x <- seq(-5, 5, 0.01)
y <- dast(x, pars = pars)
data <- rast(10^3, pars = pars)
hist(data, breaks = 50, prob = TRUE)
lines(x, y)
abline(v = 0.12, col = 2, lty = 2)

# mu
pars1 <- c(-2, 2, 0.7, 5, 5)
pars2 <- c(0, 2, 0.7, 5, 5)
pars3 <- c(2, 2, 0.7, 5, 5)
x <- seq(-5, 5, 0.01)
y1 <- dast(x, pars = pars1)
y2 <- dast(x, pars = pars2)
y3 <- dast(x, pars = pars3)
plot(x, y1, type = "l", main = expression(mu), xlab = "x", ylab = "density")
lines(x, y2, col = 2)
lines(x, y3, col = 3)
abline(v = -2, col = 1, lty = 2)
abline(v = 0, col = 2, lty = 2)
abline(v = 2, col = 3, lty = 2)
legend(x = "topleft", legend = c("mu = -2", "mu = 0", "mu = 2"),
       col = 1:3, lty = rep(1, 3))

# sigma
pars1 <- c(0, 0.5, 0.5, 5, 5)
pars2 <- c(0, 1, 0.5, 5, 5)
pars3 <- c(0, 2, 0.5, 5, 5)
x <- seq(-5, 5, 0.01)
y1 <- dast(x, pars = pars1)
y2 <- dast(x, pars = pars2)
y3 <- dast(x, pars = pars3)
plot(x, y1, type = "l", main = expression(sigma), xlab = "x", ylab = "density")
lines(x, y2, col = 2)
lines(x, y3, col = 3)
abline(v = 0, lty = 2)
legend(x = "topleft", legend = c("sigma = 0.5", "sigma = 1", "sigma = 2"),
       col = 1:3, lty = rep(1, 3))

# nu1
pars1 <- c(0, 2, 0.5, 1, 5)
pars2 <- c(0, 2, 0.5, 2, 5)
pars3 <- c(0, 2, 0.5, 5, 5)
pars4 <- c(0, 2, 0.5, 9, 5)
x <- seq(-5, 5, 0.01)
y1 <- dast(x, pars = pars1)
y2 <- dast(x, pars = pars2)
y3 <- dast(x, pars = pars3)
y4 <- dast(x, pars = pars4)
plot(x, y1, type = "l", main = expression(nu[1]), xlab = "x", ylab = "density")
lines(x, y2, lty = 2)
lines(x, y3, lty = 3)
lines(x, y4, lty = 4)
abline(v = 0, col = 4, lty = 2)
legend(x = "topleft", legend = c("nu1 = 1", "nu1 = 2", "nu1 = 5", "nu1 = 9"), lty = 1:4)

# nu2
pars1 <- c(0, 2, 0.5, 5, 1)
pars2 <- c(0, 2, 0.5, 5, 2)
pars3 <- c(0, 2, 0.5, 5, 5)
pars4 <- c(0, 2, 0.5, 5, 9)
x <- seq(-5, 5, 0.01)
y1 <- dast(x, pars = pars1)
y2 <- dast(x, pars = pars2)
y3 <- dast(x, pars = pars3)
y4 <- dast(x, pars = pars4)
plot(x, y1, type = "l", main = expression(nu[2]), xlab = "x", ylab = "density")
lines(x, y2, lty = 2)
lines(x, y3, lty = 3)
lines(x, y4, lty = 4)
abline(v = 0, col = 4, lty = 2)
legend(x = "topleft", legend = c("nu2 = 1", "nu2 = 2", "nu2 = 5", "nu2 = 9"), lty = 1:4)

# moments
pars <- c(0.12, 0.6, 0.6, 5, 6)
moment_ast(1, pars = pars)
moment_ast(2, pars = pars)
mean_ast(pars = pars)
var_ast(pars = pars)
moments_ast(pars = pars)

# fit
pars <- c(0.12, 0.6, 0.6, 6, 5)
data <- rast(1000, pars = pars)
solver_control <- list(eval.max = 10^3, iter.max = 10^3)
fit <- astfit(data, solver = 'nlminb', solver_control = solver_control)
summary(fit)
moments(fit)
plot(fit, "density")
plot(fit, "qqplot")

# simulation
data <- list(rast(10^4, 0.12, 0.6, 0.6, 3, 5))
simNlm <- my_report(data, solver = "nlminb",
                    solver_control = list(eval.max = 10^3, iter.max = 10^3), plot = FALSE)
simSol <- my_report(data, solver = "Rsolnp",
                    solver_control = list(trace = 0), plot = FALSE)
simNlo <- my_report(data, solver = "nloptr",
                    solver_control = list(algorithm = "NLOPT_LN_BOBYQA", 'maxeval' = 1.0e5, 'xtol_rel' = 1.0e-8, 'xtol_abs' = 1.0e-8), plot = FALSE)
simNlm
simSol
simNlo

# Weely returns for Large & Small cap stocks
data(retLW) # Large Weekly
data(retSW) # Small Weekly
plot.zoo(retLW)
plot.zoo(retSW)
retLW1 <- as.xts(retLW)['/2006-01-31']
retLW2 <- as.xts(retLW)['2006-01-31/']

retSW1 <- as.xts(retSW)['/2006-01-31']
retSW2 <- as.xts(retSW)['2006-01-31/']
```
```{r, stock1, eval = TRUE}
resLnlm2 <- my_report(retLW2, solver = "nlminb",
                      solver_control = list(eval.max = 10^3, iter.max = 10^3), plot = FALSE)
resLsol2 <- my_report(retLW2, solver = "Rsolnp",
                      solver_control = list(trace = 0), plot = FALSE)
resLnlo2 <- my_report(retLW2, solver = "nloptr",
                      solver_control = list(algorithm = "NLOPT_LN_BOBYQA", 'maxeval' = 1.0e5, 'xtol_rel' = 1.0e-8, 'xtol_abs' = 1.0e-8), plot = FALSE)
head(resLnlm2)
head(resLsol2)
head(resLnlo2)
resLnlm2["XOM", 1:7]
resLnlm2["XOM", 8]
resLsol2["XOM", ]

resLnlm <- my_report(retLW, solver = "nlminb",
                     solver_control = list(eval.max = 10^3, iter.max = 10^3))
resLnlmN <- my_report(retLW, solver = "nlminb",
                      solver_control = list(eval.max = 10^3, iter.max = 10^3), dist = "normal")
resLnlm1 <- my_report(retLW1, solver = "nlminb",
                      solver_control = list(eval.max = 10^3, iter.max = 10^3))
resLnlm2 <- my_report(retLW2, solver = "nlminb",
                      solver_control = list(eval.max = 10^3, iter.max = 10^3))

resSnlm <- my_report(retSW, solver = "nlminb",
                     solver_control = list(eval.max = 10^3, iter.max = 10^3))
resSnlmN <- my_report(retSW, solver = "nlminb",
                      solver_control = list(eval.max = 10^3, iter.max = 10^3), dist = "normal")
resSnlm1 <- my_report(retSW1, solver = "nlminb",
                      solver_control = list(eval.max = 10^3, iter.max = 10^3))
resSnlm2 <- my_report(retSW2, solver = "nlminb",
                      solver_control = list(eval.max = 10^3, iter.max = 10^3))
'xtol_rel' = 1.0e-8, 'xtol_abs' = 1.0e-8))


data(retSW) # Small Weekly
plot.zoo(retSW)

resSWnlm <- my_report(retSW, solver = "nlminb",
                      solver_control = list(eval.max = 10^3, iter.max = 10^3), plot = "qqplot")
resSWsol <- my_report(retSW, solver = "Rsolnp",
                      solver_control = list(trace = 0), plot = "none")
resSWnlo <- my_report(retSW, solver = "nloptr",
                      solver_control = list(algorithm = "NLOPT_LN_BOBYQA", 'maxeval' = 1.0e5, 'xtol_rel' = 1.0e-8, 'xtol_abs' = 1.0e-8), plot = "none")
