rm(list = ls())
library(nloptr)
library(st)
source("R/helper_functions.R")

seed = 22
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
resG <- nloptr::nloptr(x0,
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
# 4. check the infoMat functions
# 5. cross-sectional & real data
# 7. gatfit
# 8.



llast <- function(pars, arglist) {
  y <- arglist$data
  fixed_pars <- arglist$fixed_pars

  est_idx <- which(is.na(fixed_pars))
  fix_idx <- which(!(is.na(fixed_pars)))
  all_pars <- c()
  all_pars[est_idx] <- pars
  all_pars[fix_idx] <- fixed_pars[fix_idx]

  mu <- all_pars[1]
  sigma <- all_pars[2]
  alpha <- all_pars[3]
  nu1 <- all_pars[4]
  nu2 <- all_pars[5]
  T_ <- length(y)
  y1 <- y[y <= Re(mu)]
  y2 <- y[y > Re(mu)]

  logl <- -T_ * log(sigma) - 0.5 * (nu1 + 1) * sum(log(1 + ((y1 - mu)/(2 * alpha * sigma * K(nu1)))^2/nu1)) -
    0.5 * (nu2 + 1) * sum(log(1 + ((y2 - mu)/(2 * (1 - alpha) * sigma * K(nu2)))^2/nu2))
  -logl
}

K <- function(nu) {
  # there is a precision error here, numerical expert needed may need a c++ version
  # nu = 1e10 and greater start to be varying
  if (Re(nu) < 1000) {
    gammaz(0.5 * (nu + 1)) / sqrt(pi * nu) / gammaz(0.5 * nu)
  } else {
    # seems to be correct, need to check
    dnorm(0)
  }
}


pars = c(0, 1.5, 0.5, 2, 2)
cbind(llast_grad(pars, arglist),numDeriv::grad(llast, pars,
                                               arglist=arglist))



# newton raphson
pfunc <- function(q) {
  pgat(q, 1.5, 2, 2, 2, 2, 5) - 0.5440223
}

newton_raphson(pfunc)
x0 = 0
maxiter = 100
tol = 1e-8
h <- 1e-7
i <- 1
x1 <- x0
# p <- numeric(maxiter)
while(i <= maxiter) {
  fprime <- (pfunc(x0 + h) - pfunc(x0)) / h
  x1 <- x0 - pfunc(x0) / fprime
  # p = x1
  i <- i + 1
  if (abs(x1 - x0) < tol) {
    break
  }
  x0 <- x1
}
x1


res6 <- Rsolnp::solnp(pars = ipars$start_pars,
                      fun = llast,
                      LB= ipars$lb,
                      UB = ipars$ub,
                      arglist = arglist)



set.seed(seed)
mu <- pars[1]
sigma <- pars[2]
alpha <- pars[3]
nu1 <- pars[4]
nu2 <- pars[5]
data <- rgat(n, mu, sigma, alpha, nu1, nu2)
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
tt = Sys.time()
res3 <- nlminb(start = ipars$start_pars,
               objective = llast,
               gradient = llast_grad,
               arglist = arglist,
               control = list(eval.max = 10^3, iter.max = 10^3),
               scale = scale,
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


hist(rgat(10000, 1.5, 2, 2, 2, 2, 5), breaks = 100, freq = TRUE)
lines(seq(-100, 10, 0.1), dgat(seq(-100, 10, 0.1), 1.5, 2, 2, 2, 2, 5), type = "l")
