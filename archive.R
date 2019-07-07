## apart from astspec & information matrix for symmetric t, these are just random notes
## will clean in a future time


# # old ast-spec
# #' @title AST Specification
# #' @description Method for creating an AST distribution object prior to fitting.
# #' @param data A univariate data object, can be ... for the AST distribution to fit on.
# #' @param start_pars Numeric vector of starting parameters for the optimization algorithm.
# #' @param fixed_pars Numeric vector of parameters to be kept fixed during the optimization routine.
# #' @name astspec
# #' @examples
# #' data <- rast(1000, 0.12, 0.6, 0.7, 3, 5)
# #' spec <- astspec(data)
#
# #' @rdname astspec
# #' @export
astspec <- function(data, start_pars = c("mu" = 0, "sigma" = 1, "alpha" = 0.5, "nu1" = 1, "nu2" = 1),
                    fixed_pars = c()) {
  if (!is.numeric(data))
    stop("data must be numeric")
  if (length(start_pars) != 5)
    stop("start_pars must be a numeric of length 5")
  # if (!is.list(fixed_pars)) stop('fixed_pars must be a named list')

  check_bound(start_pars)
  bounds <- data.frame(name = c("mu", "sigma", "alpha", "nu1", "nu2"),
                       lower_bound = c(-Inf, 0, 0, 0, 0),
                       upper_bound = c(Inf, Inf, 1, Inf, Inf))
  sp_df <- data.frame(start_pars = start_pars,
                      name = names(start_pars))
  fp_df <- data.frame(fixed_pars = fixed_pars,
                      name = names(fixed_pars))
  if (length(fp_df != 0)) {
    p_df <- merge(sp_df, fp_df, by = "name", all = T)
  } else {
    sp_df$fixed_pars = NA
    p_df = sp_df
  }
  ipars <- merge(p_df, bounds, by = "name", all = T)
  rownames(ipars) <- ipars$name
  ipars$name <- NULL

  # if (length(fixed_pars) != 0) do.call(check_bound, as.list(fixed_pars))
  structure(list(data = data, ipars = ipars), class = "astspec")
}


# # Clean Session
rm(list = ls())
library(devtools)
dev_mode(on = T)
# Session Info
sessionInfo()
# remove dev packages
pkgs = installed.packages()
dev_pkgs = pkgs[, 1][pkgs[,2] == '/Users/zhiyexia/R-dev']
remove.packages(dev_pkgs)

## oo notes
# jia class
a = c(1,2,3)
is.object(a)
# in pryr
# otype for objects
# ftype for functions
# some functions do dispatch in R using UseMethod()
# while some other dipatches in C, they're called internal generics
methods(mean)
methods(class = 'ts')
# (Apart from methods defined in the base package, most S3 methods will not be visible: use getS3method() to read their
# source code.)
test <- function(x) { structure(list(), class = 'test') }
runner = test()
class(runner)
t.test <- function(x) 'Class test'
t.test(runner)
y <- 1
g <- function(x) { y <- 2
UseMethod('g') }
g.numeric <- function(x) y
g(10)
h <- function(x) { x <- 10
UseMethod('h') }
h.character <- function(x) paste('char', x)
h.numeric <- function(x)
  paste('num', x)
h('a')
h(10)
f <- function() 'this is a foo func'
f1 <- list(a = c(1,2), b = function(x) 'foo')
g <- function() 2
class(g) <- 'function'
class(f)
class(f1)
class(g)
length.function <- function(x) 'function'
length(f)
length(f1)
length(g)
library(stats4)
# From example
(mle) y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
nLL <- function(lambda) - sum(dpois(y, lambda, log = TRUE))
fit <- mle(nLL, start = list(lambda = 5), nobs = length(y))
# An S4 object
isS4(fit)
is(fit)
getClasses()
getGenerics()



# # test file # optimizer # example test # f(x, y, z) = 3 * x^3 - 2x^2 + xz - 2y^2 # x in (1,6) # y in (-3, 1) # z in (-2,
# 5) fn = function(pars) { x = pars[1] y = pars[2] z = pars[3] } # example 2 f <- function(x) 2*(x[1]-1)^2 + 5*(x[2]-3)^2 +
# 10 r <- optim(c(1, 1), f) r$convergence == 0 r$par r$value # ast problem # sigma > 0 # alpha in (0,1) # nu1 and nu2 > 0
# library(nloptr) ## Rosenbrock Banana function eval_f <- function(x) { return( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2
# ) } ## Gradient of Rosenbrock Banana function eval_grad_f <- function(x) { return( c( -400 * x[1] * (x[2] - x[1] * x[1])
# - 2 * (1 - x[1]), 200 * (x[2] - x[1] * x[1]) ) ) } x0 <- c( -1.2, 1 ) opts <- list('algorithm'='NLOPT_LD_LBFGS',
# 'xtol_rel'=1.0e-8) # solve Rosenbrock Banana function res <- nloptr( x0=x0, eval_f=eval_f, eval_grad_f=eval_grad_f,
# opts=opts) print( res ) ## Rosenbrock Banana function and gradient in one function eval_f_list <- function(x) {
# common_term <- x[2] - x[1] * x[1] return( list( 'objective' = 100 * common_term^2 + (1 - x[1])^2, 'gradient' = c( -400 *
# x[1] * common_term - 2 * (1 - x[1]), 200 * common_term) ) ) } res <- nloptr( x0=x0, eval_f=eval_f_list, opts=opts) print(
# res ) K = function(nu) { gamma(0.5*(nu+1)) / (sqrt(pi*nu)*gamma(0.5*nu)) } llast = function(pars, y) { mu = pars[1] sigma
# = pars[2] alpha = pars[3] nu1 = pars[4] nu2 = pars[5] T_ = length(y) y1 = y[y <= mu] y2 = y[y > mu] logl = T_ *
# log(sigma) + 0.5 * (nu1 + 1) * sum(log(1 + ((y1 - mu)/(2 * alpha * sigma * K(nu1)))^2/nu1)) + 0.5 * (nu2 + 1) * sum(log(1
# + ((y2 - mu)/(2 * (1-alpha) * sigma * K(nu2)))^2/nu2)) logl } y = rast(1000, 1.5, 1, 0.8, 3, 4) llast(c(1.5, 1, 0.8, 3,
# 4), y) restest = nloptr(x0 = c(0,1,0.5,1,1), eval_f = llast, lb = c(-Inf, 0, 0, 0, 0), ub = c(Inf, Inf, 1, Inf, Inf),
# opts = list('algorithm'='NLOPT_LN_COBYLA', 'maxeval' = 100000, 'xtol_rel'=1.0e-8), y = y) restest

############################################################
# with the ln_cobyla algo
# 1. the dof parameter at the side of skewness would be a problem
# 2. sigma and 2 dof would be more inaccurate
# 3. requires more time with dof smaller than 1
# 4. the est. time seems to be less when dof is larger
# 5. (1, 1.5) seems to be an exception for no.4
# 6. (1, 2) too
# 7.
############################################################


infoMat_t <- function(sigma, nu) {
  I11 <- 1/4*(trigamma((nu+1)/2) - trigamma(nu/2)) - 1/nu*(1/(nu+1) - 1/(2*(nu+3)) )
  I12 <- 1/sigma*(1/(nu+3)-1/(nu+1))
  I22 <- 2/sigma^2*nu/(nu+3)
  I13 <- I23 <- 0
  I33 <- 1/sigma^2*(1-2/(nu+3))
  infoMat <- matrix(c(I11, I12, I13,
                      I12, I22, I23,
                      I13, I23, I33),
                    nrow = 3, ncol = 3)
  rownames(infoMat) = colnames(infoMat) = c("nu", "sigma", "mu")
  infoMat
}



