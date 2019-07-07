rm(list = ls())
library(nloptr)
library(st)
source("R/helper_functions.R")

set.seed(22)
data <- rast(1000, 0.12, 0.6, 0.4, 3, 5)
ipars <- data.frame(name = c("alpha", "mu", "nu1", "nu2", "sigma"),
                    lb = c(0, -Inf, 1e-7, 1e-7, 1e-7),
                    ub = c(1, Inf, Inf, Inf, Inf),
                    start_pars = c(0.5, 0, 2, 2, 1),
                    fixed_pars = c(NA, NA, NA, NA, NA))
rownames(ipars) <- ipars$name
arglist <- list(data = data, ipars = ipars)

lb <- ipars$lb
ub <- ipars$ub
x0 <- ipars$start_pars

opts_ng <- list('algorithm' = 'NLOPT_LN_BOBYQA',
                'maxeval' = 1.0e5,
                'xtol_rel' = 1.0e-8)
resNG <- nloptr(x0,
                eval_f = llast,
                lb = lb,
                ub = ub,
                opts = opts_ng,
                arglist = arglist)
resNG
resNG$solution

res_BOBYQA <- resNG$solution
nloptr.print.options(opts.show = "algorithm")

opts_g <- list('algorithm' = 'NLOPT_LD_MMA',
             'maxeval' = 1.0e5,
             'xtol_rel' = 1.0e-8)
resG <- nloptr(x0,
               eval_f = llast,
               eval_grad_f = llast_grad,
               lb = lb,
               ub = ub,
               opts = opts_g,
               arglist = arglist)
resG
resG$solution

# rownames(ipars) <- name = c("mu", "sigma", "alpha", "nu1", "nu2")

###################
# 1. solver performance in different values
# 2. gradient based algorithms
###################

ttt1 <- function(seed, n) {
  set.seed(seed)
  data <- rast(10^n, 1.5, 2, 0.8, 0.7, 0.7)
  ipars <- data.frame(lb = c(-Inf, 0, 0, 0, 0),
                      ub = c(Inf, Inf, 1, Inf, Inf),
                      start_pars = c(0, 1, 0.5, 2, 2),
                      fixed_pars = c(NA, NA, NA, NA, NA))
  # rownames(ipars) <- name = c("mu", "sigma", "alpha", "nu1", "nu2")

  arglist <- list(data = data, fixed_pars = ipars$fixed_pars)
  lb <- ipars$lb
  ub <- ipars$ub
  x0 <- ipars$start_pars
  tt = Sys.time()
  resGL = Rsolnp::solnp(pars = x0,
                fun = llast,
                LB= lb,
                UB = ub,
                arglist = arglist)
  list(resGL$values[length(resGL$values)], resGL$pars,
  Sys.time() - tt)
}

ttt <- function(seed, n) {
  set.seed(seed)
  data <- rast(10^n, 1.5, 2, 0.8, 0.7, 0.7)
  ipars <- data.frame(lb = c(-Inf, 0, 0, 0, 0),
                      ub = c(Inf, Inf, 1, Inf, Inf),
                      start_pars = c(0, 1, 0.5, 2, 2),
                      fixed_pars = c(NA, NA, NA, NA, NA))
  # rownames(ipars) <- name = c("mu", "sigma", "alpha", "nu1", "nu2")

  arglist <- list(data = data, fixed_pars = ipars$fixed_pars)
  lb <- ipars$lb
  ub <- ipars$ub
  x0 <- ipars$start_pars
  tt = Sys.time()
  opts_ng <- list('algorithm' = 'NLOPT_LN_BOBYQA',
                  'maxeval' = 1.0e5,
                  'xtol_rel' = 1.0e-8)
  resGL <- nloptr(x0 = x0,
                  eval_f = llast,
                  lb = lb,
                  ub = ub,
                  opts = opts_ng,
                  arglist = arglist)
  list(resGL$objective, resGL$solution,
       Sys.time() - tt)
}

tt = Sys.time()
lb = c(-Inf, 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8)
opts_l <- list('algorithm' = 'NLOPT_LD_MMA',
               'maxeval' = 1.0e4, 'ftol_rel' = 1.0e-8,
               'xtol_rel' = 1.0e-8, "print_level" = 3, "check_derivatives" = TRUE)
resl <- nloptr(c(0.09, 0.6, 0.4, 3, 4.5),
               eval_f = llast,
               eval_grad_f = llast_grad,
               lb = lb,
               ub = ub,
               opts = opts_l,
               arglist = arglist)
llast_grad(resl$solution, arglist)
Sys.time() - tt


tt = Sys.time()
resGL <- nlminb(start = x0,
                objective = llast,
                gradient = llast_grad,
                arglist = arglist,
                lower = lb,
                upper = ub)
Sys.time() - tt

llast(x0, arglist)
llast_grad(x0, arglist)

# LBFGS

aa = Rsolnp::solnp(pars = x0,
              fun = llast,
              LB= lb,
              UB = ub,
              arglist = arglist)

lbfgs(call_eval = llast,
      call_grad = llast_grad,
      arglist = arglist)

optim(par = x0,
       fn = llast,
       gr = llast_grad,
      arglist = arglist,
       lower = lb,
       upper = ub,
       method = "L-BFGS-B",
       control = list(trace = 1))

# BOBYQA, MMA
nloptr.print.options(opts.show = "algorithm")

seed = 22
n = 3
set.seed(seed)
data <- rast(10^n, 1.5, 2, 0.8, 0.7, 0.7)
ipars <- data.frame(lb = c(-Inf, 0, 0, 0, 0),
                    ub = c(Inf, Inf, 1, Inf, Inf),
                    start_pars = c(0, 1, 0.5, 2, 2),
                    fixed_pars = c(NA, NA, NA, NA, NA))
# rownames(ipars) <- name = c("mu", "sigma", "alpha", "nu1", "nu2")
arglist <- list(data = data, fixed_pars = ipars$fixed_pars)
lb <- ipars$lb
ub <- ipars$ub
x0 <- ipars$start_pars
tt = Sys.time()
resGL = Rsolnp::solnp(pars = x0,
                      fun = llast,
                      LB= lb,
                      UB = ub,
                      arglist = arglist)
list(resGL$values[length(resGL$values)], resGL$pars,
     Sys.time() - tt)

# 1. reproducible
# 2. authorized domain
# 3. mean, var, ...
