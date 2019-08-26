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

################################################
# this belongs to astfit
################################################



# with the grid search commented, it is now just a wrapper
# # grid search 1
# if ("nu1" %in% names(fixed_pars)) { nu1vec <- fixed_pars["nu1"] } else { nu1vec <- seq(2, 20, by = 4) }
# if ("nu2" %in% names(fixed_pars)) { nu2vec <- fixed_pars["nu2"] } else { nu2vec <- seq(2, 20, by = 4) }
# grid <- array(c(rep(nu1vec, length(nu2vec)), rep(nu2vec, length(nu1vec))), c(length(nu1vec), length(nu2vec), 2))
# if (dim(grid)[1] == dim(grid)[2]) { grid[,,2] = t(grid[,,2]) }
# valueGrid <- apply(grid, 1:2, objective_value, data, start_pars, fixed_pars, solver, solver_control)
# idx <- as.numeric(which(valueGrid == min(valueGrid), arr.ind = T))
# idx <- c(grid[idx[1], idx[2], 1], grid[idx[1], idx[2], 2])
#
# # grid search 2
# if (!("nu1" %in% names(fixed_pars))) { nu1vec <- seq(idx[1]-2, idx[1]+2, by=1) }
# if (!("nu2" %in% names(fixed_pars))) { nu2vec <- seq(idx[2]-2, idx[2]+2, by=1) }
# grid <- array(c(rep(nu1vec, length(nu2vec)), rep(nu2vec, length(nu1vec))), c(length(nu1vec), length(nu2vec), 2))
# if (dim(grid)[1] == dim(grid)[2]) { grid[,,2] = t(grid[,,2]) }
# valueGrid <- apply(grid, 1:2, objective_value, data, start_pars, fixed_pars, solver, solver_control)
# idx <- as.numeric(which(valueGrid == min(valueGrid), arr.ind = T))
# idx <- c(grid[idx[1], idx[2], 1], grid[idx[1], idx[2], 2])
#
# # fit with result of grid search 2
# fp_tmp <- fixed_pars
# fp_tmp["nu1"] <- idx[1]
# fp_tmp["nu2"] <- idx[2]
# fit_tmp <- astfit_local(data, start_pars, fp_tmp, solver, solver_control)
#
# # final fit
# start_pars <- c(fit_tmp$sol_res$solution, c(idx[1], idx[2]))
# names(start_pars) <- c("mu", "sigma", "alpha", "nu1", "nu2")

#### this too ####
objective_value <- function(nus, data, start_pars, fixed_pars, solver, solver_control) {
  if (nus[1] == 0 || nus[2] == 0) {
    obj <- 10^5
  } else {
    fixed_pars["nu1"] = nus[1]
    fixed_pars["nu2"] = nus[2]
    obj <- astfit_local(data, start_pars, fixed_pars, solver, solver_control)$sol_res$objective
  }
  obj
}

llast <- function(y, mu, sigma, alpha, nu1, nu2) {
  T_ <- length(y)
  y1 <- y[y <= mu]
  y2 <- y[y > mu]

  logl <- -T_ * log(sigma) - 0.5 * (nu1 + 1) * sum(log(1 + ((y1 - mu)/(2 * alpha * sigma * K(nu1)))^2/nu1)) -
    0.5 * (nu2 + 1) * sum(log(1 + ((y2 - mu)/(2 * (1 - alpha) * sigma * K(nu2)))^2/nu2))
  -logl
}


mean_gat <- function(mu, phi, alpha, r, c, nu) {
  A <- nu/( alpha*(1+r^2) )
  B <- A * r^2
  delta <- A / nu * r
  mu + phi * (beta(A + delta, B - delta)/c - c*beta(A - delta, B + delta)) /
    (2*beta(A, B))
}

var_gat <- function(mu, phi, alpha, r, c, nu) {
  A <- nu/( alpha*(1+r^2) )
  B <- A * r^2
  delta <- A / nu * r
  phi^2 / (4*beta(A, B)) * (c^(-2)*beta(A + 2*delta, B - 2*delta) + c^2*beta(A - 2*delta, B + 2* delta)) -
    phi^2/2
}

# var_gat = moment_central_gat but not var(data)
vg <- function(mu, phi, alpha, r, c, nu) {
  m <- mean_gat(mu, phi, alpha, r, c, nu)
  integrand <- function(x) {
    (x - m)^2 * dgat(x, mu, phi, alpha, r, c, nu)
  }
  safeIntegrate(integrand, -Inf, Inf)$value
}

sg <- function(mu, phi, alpha, r, c, nu) {
  m <- mean_gat(mu, phi, alpha, r, c, nu)
  sd <- sqrt(vg(mu, phi, alpha, r, c, nu))
  integrand <- function(x) {
    ((x - m)/sd)^3 * dgat(x, mu, phi, alpha, r, c, nu)
  }
  safeIntegrate(integrand, -Inf, Inf)$value
}

#' @rdname gat-methods
#' @export
fitted.gat <- function(fit) {
  fit$fitted_pars
}

#' @rdname gat-methods
#' @export
se.gat <- function(fit) {
  fit$standard_errors
}

#' @rdname gat-methods
#' @export
objective.gat <- function(fit) {
  fit$objective
}


#' @export
se <- function(x, ...) {
  UseMethod("se", x)
}

#' @export
objective <- function(x, ...) {
  UseMethod("objective", x)
}

#' @export
fitted <- function(x, ...) {
  UseMethod("fitted", x)
}


# has no documentation developed yet
#' @export
surfacePlot <- function(n, pars, plotPars, ...) {
  mu <- pars[1]
  sigma <- pars[2]
  alpha <- pars[3]
  nu1 <- pars[4]
  nu2 <- pars[5]
  data <- rast(n, mu, sigma, alpha, nu1, nu2)

  xName <- plotPars[1]
  yName <- plotPars[2]
  xVec <- parVec(pars[xName], xName)
  yVec <- parVec(pars[yName], yName)
  xLen <- length(xVec)
  yLen <- length(yVec)
  xMat <- matrix(rep(xVec, yLen), xLen, yLen)
  yMat <- matrix(rep(yVec, xLen), xLen, yLen, byrow = TRUE)
  parGrid <- array(c(xMat, yMat), c(xLen, yLen, 2))

  start_pars <- c()
  fixed_pars <- c()
  solver <- "Rsolnp"
  solver_control <- list(trace = 0)
  valGrid <- apply(parGrid, 1:2, obj_surface, data, start_pars, fixed_pars, solver, solver_control, xName, yName)
  rownames(valGrid) <- xVec
  colnames(valGrid) <- yVec
  persp(xVec, yVec, valGrid, xlab = xName, ylab = yName, ...)
  return(list(xVec, yVec, valGrid))
}

report <- function(fit) {
  report <- c(round(c(fit$fitted_pars, fit$objective, fit$time_elapsed) , 8), fit$message)
  names(report)[c(6, 7, 8)] <- c("objective", "time", "message")
  report
}

# has no documentation developed yet
#' @export
my_report <- function(data, solver, solver_control, plots = "none", dist = "ast") {
  fitList <- lapply(data, astfit, solver = solver, solver_control = solver_control)

  for (i in 1:length(fitList)) {
    fitList[[i]]$name <- names(fitList)[i]
  }

  if (plots == "both") {
    par(mfrow = c(4, 5))
    lapply(fitList, function(fit) {plot(fit, type = "density", main = fit$name) })
    par(mfrow = c(4, 5))
    lapply(fitList, function(fit) {plot(fit, type = "qqplot", main = fit$name, dist = dist) })
    #}
  } else if (plots == "density") {
    par(mfrow = c(4, 5))
    lapply(fitList, function(fit) {plot(fit, type = "density", main = fit$name) })
  } else if (plots == "qqplot") {
    par(mfrow = c(4, 5))
    lapply(fitList, function(fit) {plot(fit, type = "qqplot", main = fit$name, dist = dist) })
  }
  par(mfrow = c(1, 1))

  res <- t(sapply(fitList, report))
  res <- as.data.frame(res)
  res[,-length(res)] <- lapply(res[,-length(res)], function(x) round(as.numeric(as.character(x)), 6))

  res
}

newton_raphson <- function(f, x0 = 0.5, maxiter = 1e2, tol = 1e-8) {
  h <- 1e-8
  i <- 1
  x1 <- x0
  x <- x0
  # p <- numeric(maxiter)
  while(i <= maxiter) {
    fprime <- (f(x0 + h) - f(x0)) / h
    x1 <- (x0 - (f(x0) / fprime))
    # p = x1
    i <- i + 1
    if (abs(x1 - x0) < tol) {
      break
    }
    x0 <- x1
  }
  x1
}

# used only for surfaceplot in plot-methods
obj_surface <- function(pars, data, start_pars, fixed_pars, solver, solver_control, xName, yName, symmetric = FALSE) {
  fixed_pars[xName] <- pars[1]
  fixed_pars[yName] <- pars[2]
  #astfit_local(data, start_pars, fixed_pars, solver, solver_control, symmetric)
  return(astfit_local(data, start_pars, fixed_pars, solver, solver_control, symmetric)$objective)
}

safeIntegrate <- function(f, lower, upper, ..., subdivisions = 1000L,
                          rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
                          stop.on.error = TRUE, keep.xy = FALSE, aux = NULL) {
  res <- integrate(f, lower, upper, ..., subdivisions = subdivisions,
                   rel.tol = rel.tol, abs.tol = abs.tol,
                   stop.on.error = stop.on.error, keep.xy = keep.xy, aux = aux)
  # the safeIntegrate in HyperbolicDist has different rules
  # don't seem to be necessary
  if (lower == upper) {
    res$value <- 0
    res$abs.error <- 0
  }
  res
}

# the following 2 are required for surface plot, for exploraton reason
parVec <- function(x, xName) {
  eps <- 1.0e-8
  if (xName == "mu") {
    seq(x - 0.5, x + 0.5, 0.1)
  } else if (xName == "alpha") {
    seq(0.1, 0.9, 0.1)
  } else {
    seq(max(0.1, x - 0.5), x + 0.5, 0.1)
  }
}
