rm(list = ls())
library(nloptr)
library(st)
source("R/helper-functions.R")

seed = 22
set.seed(seed)
n = 4
data <- rast(10^n, 0.12, 0.6, 0.4, 3, 5)
ipars <- data.frame(lb = c(-Inf, 0, 0, 0, 0),
                    ub = c(Inf, Inf, 1, Inf, Inf),
                    start_pars = c(0, 1, 0.5, 2, 2),
                    fixed_pars = c(NA, NA, NA, NA, NA))
# rownames(ipars) <- name = c("mu", "sigma", "alpha", "nu1", "nu2")
pars = c(0.12, 0.6, 0.4, 3, 5)
data = rast(10^4, pars)
solver_control <- list(eval.max = 10^4, iter.max = 10^4)
fit <- astfit(data, solver = 'nlminb', solver_control = solver_control)

arglist <- list(data = data, fixed_pars = ipars$fixed_pars)
llast(pars, arglist)
llast_grad(pars, arglist)

cbind(llast_grad(pars, arglist),numDeriv::grad(llast, pars,
                                               arglist=arglist))


lb <- ipars$lb
ub <- ipars$ub
x0 <- ipars$start_pars

nloptr.print.options(opts.show = "algorithm")

# 1. reproducible
# 2. authorized domain
# 3. mean, var, ...
# 4. check the infoMat functions
# 5. cross-sectional & real data
# 7. gatfit

hist(rgat(10000, 1.5, 2, 2, 2, 2, 5), breaks = 100, freq = TRUE)
lines(seq(-100, 10, 0.1), dgat(seq(-100, 10, 0.1), 1.5, 2, 2, 2, 2, 5), type = "l")




x1 <- seq(-3, 3, 0.01)

 <-
# mu
pars11 <- c(-2, 2, 0.7, 5, 5)
pars12 <- c(0, 2, 0.7, 5, 5)
pars13 <- c(2, 2, 0.7, 5, 5)
y11 <- dast(x2, pars = pars11)
y12 <- dast(x2, pars = pars12)
y13 <- dast(x2, pars = pars13)
plot(x2, y11, type = "l", main = expression(mu), xlab = "x", ylab = "density")
lines(x2, y12, col = 2)
lines(x2, y13, col = 3)
abline(v = -2, col = 1, lty = 2)
abline(v = 0, col = 2, lty = 2)
abline(v = 2, col = 3, lty = 2)
legend(x = "topleft", legend = c("mu = -2", "mu = 0", "mu = 2"), col = 1:3, lty = rep(1, 3))

# sigma
pars21 <- c(0, 0.5, 0.5, 5, 5)
pars22 <- c(0, 1, 0.5, 5, 5)
pars23 <- c(0, 2, 0.5, 5, 5)
y21 <- dast(x1, pars = pars21)
y22 <- dast(x1, pars = pars22)
y23 <- dast(x1, pars = pars23)
plot(x1, y21, type = "l", main = expression(sigma), xlab = "x", ylab = "density")
lines(x1, y22, col = 2)
lines(x1, y23, col = 3)
abline(v = 0, lty = 2)
legend(x = "topleft", legend = c("sigma = 0.5", "sigma = 1", "sigma = 2"), col = 1:3, lty = rep(1, 3))

# alpha
pars31 <- c(0, 2, 0.3, 5, 5)
pars32 <- c(0, 2, 0.5, 5, 5)
pars33 <- c(0, 2, 0.7, 5, 5)
y31 <- dast(x1, pars = pars31)
y32 <- dast(x1, pars = pars32)
y33 <- dast(x1, pars = pars33)
plot(x1, y31, type = "l", main = expression(alpha), xlab = "x", ylab = "density")
lines(x1, y32, lty = 2)
lines(x1, y33, lty = 3)
abline(v = 0, col = 4, lty = 2)
legend(x = "topleft", legend = c("alpha = 0.3", "alpha = 0.5", "alpha = 0.7"), lty = 1:3)

# nu1
pars41 <- c(0, 2, 0.5, 1, 5)
pars42 <- c(0, 2, 0.5, 2, 5)
pars43 <- c(0, 2, 0.5, 5, 5)
pars44 <- c(0, 2, 0.5, 9, 5)
y41 <- dast(x1, pars = pars41)
y42 <- dast(x1, pars = pars42)
y43 <- dast(x1, pars = pars43)
y44 <- dast(x1, pars = pars44)
plot(x1, y41, type = "l", main = expression(nu[1]), xlab = "x", ylab = "density")
lines(x1, y42, lty = 2)
lines(x1, y43, lty = 3)
lines(x1, y44, lty = 4)
abline(v = 0, col = 4, lty = 2)
legend(x = "topleft", legend = c("nu1 = 1", "nu1 = 2", "nu1 = 5", "nu1 = 9"), lty = 1:4)

# nu2
pars51 <- c(0, 2, 0.5, 5, 1)
pars52 <- c(0, 2, 0.5, 5, 2)
pars53 <- c(0, 2, 0.5, 5, 5)
pars54 <- c(0, 2, 0.5, 5, 9)
y51 <- dast(x1, pars = pars51)
y52 <- dast(x1, pars = pars52)
y53 <- dast(x1, pars = pars53)
y54 <- dast(x1, pars = pars54)
plot(x1, y51, type = "l", main = expression(nu[2]), xlab = "x", ylab = "density")
lines(x1, y52, lty = 2)
lines(x1, y53, lty = 3)
lines(x1, y54, lty = 4)
abline(v = 0, col = 4, lty = 2)
legend(x = "topleft", legend = c("nu2 = 1", "nu2 = 2", "nu2 = 5", "nu2 = 9"), lty = 1:4)

# fit
pars <- c(0.12, 0.6, 0.6, 3, 5)
data <- rast(1000, pars = pars)
solver_control <- list(eval.max = 10^3, iter.max = 10^3)
fit <- astfit(data, solver = 'nlminb', solver_control = solver_control)
summary(fit)
moments(fit)
plot(fit)

# moments
pars <- c(0.12, 0.6, 0.6, 3, 5)
moment_ast(1, pars = pars)
moment_ast(2, pars = pars)
mean_ast(pars = pars)
var_ast(pars = pars)
skew_ast(pars = pars)
kurt_ast(pars = pars)
moments_ast(pars = pars)

##### gat
x1 <- seq(-3, 3, 0.01)
x2 <- seq(-5, 5, 0.01)
x3 <- seq(-10, 10, 0.01)

pars0 <- c(1.5, 2, 0.7, 3, 5)
y0 <- dast(x2, pars = pars0)
plot(x2, y0, type = "l", xlab = "x", ylab = "y")
abline(v = 1.5, col = 4, lty = 2)


# mu
pars11 <- c(-2, 2, 0.7, 5, 5)
pars12 <- c(0, 2, 0.7, 5, 5)
pars13 <- c(2, 2, 0.7, 5, 5)
y11 <- dast(x2, pars = pars11)
y12 <- dast(x2, pars = pars12)
y13 <- dast(x2, pars = pars13)
plot(x2, y11, type = "l", main = expression(mu), xlab = "x", ylab = "density")
lines(x2, y12, col = 2)
lines(x2, y13, col = 3)
abline(v = -2, col = 1, lty = 2)
abline(v = 0, col = 2, lty = 2)
abline(v = 2, col = 3, lty = 2)
legend(x = "topleft", legend = c("mu = -2", "mu = 0", "mu = 2"), col = 1:3, lty = rep(1, 3))

# sigma
pars21 <- c(0, 1, 2, 2, 2, 5)
pars22 <- c(0, 2, 2, 2, 2, 5)
pars23 <- c(0, 3, 2, 2, 2, 5)
y21 <- dgat(x1, pars = pars21)
y22 <- dgat(x1, pars = pars22)
y23 <- dgat(x1, pars = pars23)
plot(x1, y21, type = "l", main = expression(phi), xlab = "x", ylab = "density")
lines(x1, y22, col = 2)
lines(x1, y23, col = 3)
abline(v = 0, lty = 2)
legend(x = "topleft", legend = c("sigma = 0.5", "sigma = 1", "sigma = 2"), col = 1:3, lty = rep(1, 3))

# alpha
pars31 <- c(0, 2, 0.3, 5, 5)
pars32 <- c(0, 2, 0.5, 5, 5)
pars33 <- c(0, 2, 0.7, 5, 5)
y31 <- dast(x1, pars = pars31)
y32 <- dast(x1, pars = pars32)
y33 <- dast(x1, pars = pars33)
plot(x1, y31, type = "l", main = expression(alpha), xlab = "x", ylab = "density")
lines(x1, y32, lty = 2)
lines(x1, y33, lty = 3)
abline(v = 0, col = 4, lty = 2)
legend(x = "topleft", legend = c("alpha = 0.3", "alpha = 0.5", "alpha = 0.7"), lty = 1:3)

# nu1
pars41 <- c(0, 2, 0.5, 1, 5)
pars42 <- c(0, 2, 0.5, 2, 5)
pars43 <- c(0, 2, 0.5, 5, 5)
pars44 <- c(0, 2, 0.5, 9, 5)
y41 <- dast(x1, pars = pars41)
y42 <- dast(x1, pars = pars42)
y43 <- dast(x1, pars = pars43)
y44 <- dast(x1, pars = pars44)
plot(x1, y41, type = "l", main = expression(nu[1]), xlab = "x", ylab = "density")
lines(x1, y42, lty = 2)
lines(x1, y43, lty = 3)
lines(x1, y44, lty = 4)
abline(v = 0, col = 4, lty = 2)
legend(x = "topleft", legend = c("nu1 = 1", "nu1 = 2", "nu1 = 5", "nu1 = 9"), lty = 1:4)

# nu2
pars51 <- c(0, 2, 0.5, 5, 1)
pars52 <- c(0, 2, 0.5, 5, 2)
pars53 <- c(0, 2, 0.5, 5, 5)
pars54 <- c(0, 2, 0.5, 5, 9)
y51 <- dast(x1, pars = pars51)
y52 <- dast(x1, pars = pars52)
y53 <- dast(x1, pars = pars53)
y54 <- dast(x1, pars = pars54)
plot(x1, y51, type = "l", main = expression(nu[2]), xlab = "x", ylab = "density")
lines(x1, y52, lty = 2)
lines(x1, y53, lty = 3)
lines(x1, y54, lty = 4)
abline(v = 0, col = 4, lty = 2)
legend(x = "topleft", legend = c("nu2 = 1", "nu2 = 2", "nu2 = 5", "nu2 = 9"), lty = 1:4)





n <- 10^4
pars <- c(mu = 0.12, sigma = 0.6, alpha = 0.6, nu1 = 3, nu2 = 3)
par(mfrow = c(2, 3))
ms <- surfacePlot(n, pars, c("mu", "sigma"), theta = 100, col = 3, shade = 0.75)
ma <- surfacePlot(n, pars, c("mu", "alpha"), theta = 100, col = 3, shade = 0.75)
mn1 <- surfacePlot(n, pars, c("mu", "nu1"), theta = 0, col = 3, shade = 0.75)
mn2 <- surfacePlot(n, pars, c("mu", "nu2"), theta = 0, col = 3, shade = 0.75)
sa <- surfacePlot(n, pars, c("sigma", "alpha"), theta = 40, col = 3, shade = 0.75)
sn1 <- surfacePlot(n, pars, c("sigma", "nu1"), theta = 40, col = 3, shade = 0.75)
par(mfrow = c(2, 2))
sn2 <- surfacePlot(n, pars, c("sigma", "nu2"), theta = 40, col = 3, shade = 0.75)
an1 <- surfacePlot(n, pars, c("alpha", "nu1"), theta = 140, col = 3, shade = 0.75)
an2 <- surfacePlot(n, pars, c("alpha", "nu2"), theta = 40, col = 3, shade = 0.75)
n1n2 <- surfacePlot(n, pars, c("nu1", "nu2"), theta = 70, col = 3, shade = 0.75)
