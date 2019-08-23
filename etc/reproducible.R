library(st)
source("R/helper-functions.R") # for llast & llast_grad
source("R/astfit.R") # for llast & llast_grad

# 1. BOBYQA
# 2. Rsolnp
# 3. LBFGS
# 4. nlminb

report <- function(pars, n, seed, flag = 0) {
  set.seed(seed)
  mu <- pars[1]
  sigma <- pars[2]
  alpha <- pars[3]
  nu1 <- pars[4]
  nu2 <- pars[5]
  data <- rast(n, mu, sigma, alpha, nu1, nu2)
  # 1.0e-8 instead of 0 to prevent division by 0, not very helpful though
  ipars <- data.frame(lb = c(-Inf, 1e-08, 1e-08, 1e-08, 1e-08),
                      ub = c(Inf, Inf, 1, Inf, Inf),
                      start_pars = c(0, 1, 0.5, 2, 2),
                      fixed_pars = c(NA, NA, NA, NA, NA),
                      names = c("mu", "sigma", "alpha", "nu1", "nu2"),
                      order = 1:5)
  arglist <- list(data = data, fixed_pars = ipars$fixed_pars)
  scale <- pars
  scale[1] <- ifelse(mu == 0, c(1, sigma, alpha, nu1, nu2), pars)

  if (flag == 2) {
    tt = Sys.time()
    res1 <- optim(par = ipars$start_pars,
                  fn = llast,
                  gr = llast_grad,
                  arglist = arglist,
                  method = "L-BFGS-B",
                  lower = ipars$lb,
                  upper = ipars$ub,
                  control = list(parscale = scale))
    t1 = Sys.time() - tt
    tt = Sys.time()
    res2 <- optim(par = ipars$start_pars,
                  fn = llast,
                  gr = llast_grad,
                  arglist = arglist,
                  method = "L-BFGS-B",
                  lower = ipars$lb,
                  upper = ipars$ub,
                  control = list(trace = 1))
    t2 = Sys.time() - tt
  } else if (flag == 1) {
    tt = Sys.time()
    res2 <- optim(par = ipars$start_pars,
                  fn = llast,
                  gr = llast_grad,
                  arglist = arglist,
                  method = "L-BFGS-B",
                  lower = ipars$lb,
                  upper = ipars$ub,
                  control = list(trace = 1))
    t2 = Sys.time() - tt
  }
  # tt = Sys.time()
  # res3 <- nlminb(start = ipars$start_pars,
  #                objective = llast,
  #                gradient = llast_grad,
  #                arglist = arglist,
  #                control = list(eval.max = 10^3, iter.max = 10^3),
  #                scale = scale,
  #                lower = ipars$lb,
  #                upper = ipars$ub)
  # t3 = Sys.time() - tt
  # tt = Sys.time()
  # res4 <- nlminb(start = ipars$start_pars,
  #                objective = llast,
  #                gradient = llast_grad,
  #                arglist = arglist,
  #                control = list(eval.max = 10^3, iter.max = 10^3),
  #                lower = ipars$lb,
  #                upper = ipars$ub)
  # t4 = Sys.time() - tt
  tt = Sys.time()
  res5 <- nloptr::nloptr(x0 = ipars$start_pars,
                         eval_f = llast,
                         lb = ipars$lb,
                         ub = ipars$ub,
                         opts = list('algorithm' = 'NLOPT_LN_BOBYQA', 'maxeval' = 1.0e5, 'xtol_rel' = 1.0e-8),
                         arglist = arglist)
  t5 = Sys.time() - tt
  res6 <- Rsolnp::solnp(pars = ipars$start_pars,
                        fun = llast,
                        LB= ipars$lb,
                        UB = ipars$ub,
                        arglist = arglist)
  if (flag == 2) {
    res <- as.data.frame(rbind(res1$par, res2$par, res3$par, res4$par, res5$solution, res6$pars))
    rownames(res) <- c("L-BFGS-B-scale", "L-BFGS-B", "nlminb-scale", "nlminb", "nloptr-bobyqa", "solnp")
    colnames(res) <- c("mu", "sigma", "alpha", "nu1", "nu2")
    res$objective <- c(res1$value, res2$value, res3$objective, res4$objective, res5$objective, res6$values[length(res6$values)])
    res$time <- c(t1, t2, t3, t4, t5, res6$elapsed)
    res$message <- c(res1$message, res2$message, res3$message, res4$message, res5$message, res6$convergence)
  } else if (flag == 1) {
    res <- as.data.frame(rbind(res2$par, res3$par, res4$par, res5$solution, res6$pars))
    rownames(res) <- c("L-BFGS-B", "nlminb-scale", "nlminb", "nloptr-bobyqa", "solnp")
    colnames(res) <- c("mu", "sigma", "alpha", "nu1", "nu2")
    res$objective <- c(res2$value, res3$objective, res4$objective, res5$objective, res6$values[length(res6$values)])
    res$time <- c(t2, t3, t4, t5, res6$elapsed)
    res$message <- c(res2$message, res3$message, res4$message, res5$message, res6$convergence)
  } else {
    res <- as.data.frame(rbind(res5$solution, res6$pars))
    rownames(res) <- c("nloptr-bobyqa", "solnp")
    colnames(res) <- c("mu", "sigma", "alpha", "nu1", "nu2")
    res$objective <- c(res5$objective, res6$values[length(res6$values)])
    res$time <- c(t5, res6$elapsed)
    res$message <- c(res5$message, res6$convergence)
  }
  res
}

seed = 22
n = 10^3
pars = c(mu = 0.12, sigma = 0.6, alpha = 0.3, nu1 = 3, nu2 = 5)
report(pars, n, seed)

# scale of mu # and sigma
# put it in an equation
seed = 22
n = 10^4
#pars1 = c(150, 2, 0.9, 2, 4)
# pars1 = c(1.5, 2, 0., 2, 4)
report(pars, n, seed, flag = 0)


## function - taylor series
## lgamma gamma gammaz



report(pars1, n, seed, flag = 1)
report(pars1, n, seed, flag = 2)
# data = rast(10^3, 1, 20, 0.99, 2, 4)
# arglist = list(data = data, fixed_pars = rep(NA, 5))
# llast(pars = pars1, arglist = arglist)

pars2 = c(1.5, 2, 0.8, 0.7, 0.7)
report(pars1, n, seed)

# proof that gradient is most likely correct
data = rast(1000, 0.12, 0.6, 0.3, 3, 5)
pars = c(0, 1.5, 0.5, 2, 2)
arglist = list(data = data, fixed_pars = rep(NA, 5))
cbind(llast_grad(pars, arglist),numDeriv::grad(llast, pars,
                                               arglist=arglist))


llast(c(0.12, 0.6, 0.4, 3 + complex(1), 5), arglist)
# surface plot
surf = surfacePlot(1000, pars, c("nu1", "nu2"), theta = 100, col = 3, shade = 0.75)
