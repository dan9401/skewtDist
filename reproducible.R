library(st)
source("R/helper-functions.R") #
source("R/astfit.R") # for llast & llast_grad

# 1. BOBYQA
# 2. Rsolnp
# 3. LBFGS
# 4. nlminb

report <- function(pars, n, seed) {
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

  # tt = Sys.time()
  # res1 <- optim(par = ipars$start_pars,
  #              fn = llast,
  #              gr = llast_grad,
  #              arglist = arglist,
  #              method = "L-BFGS-B",
  #              lower = ipars$lb,
  #              upper = ipars$ub,
  #              control = list(parscale = pars))
  # t1 = Sys.time() - tt
  # tt = Sys.time()
  # res2 <- optim(par = ipars$start_pars,
  #               fn = llast,
  #               gr = llast_grad,
  #               arglist = arglist,
  #               method = "L-BFGS-B",
  #               lower = ipars$lb,
  #               upper = ipars$ub,
  #               control = list(trace = 1))
  # t2 = Sys.time() - tt
  tt = Sys.time()
  res3 <- nlminb(start = ipars$start_pars,
                 objective = llast,
                 gradient = llast_grad,
                 arglist = arglist,
                 control = list(eval.max = 10^3, iter.max = 10^3),
                 scale = pars,
                 lower = ipars$lb,
                 upper = ipars$ub)
  t3 = Sys.time() - tt
  tt = Sys.time()
  res4 <- nlminb(start = ipars$start_pars,
                 objective = llast,
                 gradient = llast_grad,
                 arglist = arglist,
                 control = list(eval.max = 10^3, iter.max = 10^3),
                 lower = ipars$lb,
                 upper = ipars$ub)
  t4 = Sys.time() - tt
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
  # res <- as.data.frame(rbind(res1$par, res2$par, res3$par, res4$par, res5$solution, res6$pars))
  # rownames(res) <- c("L-BFGS-B-scale", "L-BFGS-B", "nlminb-scale", "nlminb", "nloptr-bobyqa", "solnp")
  # colnames(res) <- c("mu", "sigma", "alpha", "nu1", "nu2")
  # res$objective <- c(res1$value, res2$value, res3$objective, res4$objective, res5$objective, res6$values[length(res6$values)])
  # res$time <- c(t1, t2, t3, t4, t5, res6$elapsed)
  # res$message <- c(res1$message, res2$message, res3$message, res4$message, res5$message, res6$convergence)

  res <- as.data.frame(rbind(res3$par, res4$par, res5$solution, res6$pars))
  rownames(res) <- c("nlminb-scale", "nlminb", "nloptr-bobyqa", "solnp")
  colnames(res) <- c("mu", "sigma", "alpha", "nu1", "nu2")
  res$objective <- c(res3$objective, res4$objective, res5$objective, res6$values[length(res6$values)])
  res$time <- c(t3, t4, t5, res6$elapsed)
  res$message <- c(res3$message, res4$message, res5$message, res6$convergence)
  res
}

seed = 22
n = 10^3
pars = c(mu = 0.12, sigma = 0.6, alpha = 0.3, nu1 = 3, nu2 = 5)
report(pars, n, seed)

pars1 = c(1.5, 2, 0.8, 2, 4)
report(pars1, n, seed)

pars2 = c(1.5, 2, 0.8, 2, 4)
report(pars1, n, seed)

# proof that gradient is most likely correct
data = rast(10000, 0.12, 0.6, 0.3, 3, 5)
arglist = list(data = data, fixed_pars = rep(NA, 5))
llast_grad(c(0, 1, 0.5, 2, 2), arglist)
c((llast(c(0.001, 1, 0.5, 2, 2), arglist) - llast(c(-0.001, 1, 0.5, 2, 2), arglist))/0.002,
  (llast(c(0, 1.001, 0.5, 2, 2), arglist) - llast(c(0, 0.999, 0.5, 2, 2), arglist))/0.002,
  (llast(c(0, 1, 0.501, 2, 2), arglist) - llast(c(0, 1, 0.499, 2, 2), arglist))/0.002,
  (llast(c(0, 1, 0.5, 2.001, 2), arglist) - llast(c(0, 1, 0.5, 1.999, 2), arglist))/0.002,
  (llast(c(0, 1, 0.5, 2, 2.001), arglist) - llast(c(0, 1, 0.5, 2, 1.999), arglist))/0.002)
