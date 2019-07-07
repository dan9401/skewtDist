# may keep separate files for gat & ast distribuitons

# these 4 are quality of life functions, to keep the formulas simple, the list may be extended
# as the current hessian is still very complicated for ast, not mentioning the unimplemented gat
K <- function(nu) {
  # there is a precision error here, numerical expert needed may need a c++ version
  # nu = 1e10 and greater start to be varying
  if (nu < 1000) {
    exp(lgamma(0.5 * (nu + 1)) - log(sqrt(pi * nu)) - lgamma(0.5 * nu))
  } else {
    # seems to be correct, need to check
    dnorm(0)
  }
}

D <-  function(nu) {
  D <- digamma((nu + 1) / 2) - digamma(nu / 2)
  D
}

Dprime <- function(nu) {
  Dp <- trigamma((nu+1)/2) - trigamma(nu/2)
}

L <- function(pars, y) {
  mu <- pars[1]
  alpha <- pars[3]
  nu1 <- pars[4]
  L <- 1 + 1 / nu1 * ( (y - mu) / (2 * alpha * K(nu1)) )^2
  L
}

R <- function(pars, y) {
  mu <- pars[1]
  alpha <- pars[3]
  nu2 <- pars[5]
  R <- 1 + 1 / nu2 * ( (y - mu) / (2 * (1 - alpha) * K(nu2)) )^2
  R
}

# the following 2 are required for surface plot, for exploraton reason
parVec <- function(x, xName) {
  eps <- 1.0e-8
  if (xName == "mu") {
    parVal(x, -Inf, Inf)
  } else if (xName == "alpha") {
    parVal(x, eps, 1)
  } else {
    parVal(x, eps, Inf)
  }
}

parVal <- function(x, lower, upper) {
  if (abs(x) < 1) {
    seq(max(x - 0.5, lower), min(x + 0.5, upper), length.out = 11)
  } else if (abs(x) < 4) {
    seq(max(x - 1, lower), min(x + 1, upper), length.out = 11)
  } else {
    seq(max(x - 2, lower), min(x + 2, upper), length.out = 11)
  }
}

# this may be discarded or put into another form, not sure at the time
# but most likely needed
check_bound <- function(pars) {
  mu <- pars["mu"]
  sigma <- pars["sigma"]
  alpha <- pars["alpha"]
  nu1 <- pars["nu1"]
  nu2 <- pars["nu2"]

  if (!is.numeric(mu))
    stop("mu must be numeric")
  if (!is.numeric(sigma))
    stop("sigma must be numeric")
  if (!is.numeric(alpha))
    stop("alpha must be numeric")
  if (!is.numeric(nu1))
    stop("nu1 must be numeric")
  if (!is.numeric(nu2))
    stop("nu2 must be numeric")
  if (sigma <= 0)
    stop("sigma must be greater than 0")
  if (nu1 <= 0)
    stop("nu1 must be greater than 0")
  if (nu2 <= 0)
    stop("nu2 must be greater than 0")
  if (alpha <= 0 || alpha >= 1)
    stop("alpha must be between 0 and 1")
}

# used only for surfaceplot in plot-methods
obj_surface <- function(pars, data, start_pars, fixed_pars, solver, solver_control, xName, yName) {
  fixed_pars[xName] <- pars[1]
  fixed_pars[yName] <- pars[2]
  return(astfit_local(data, start_pars, fixed_pars, solver, solver_control)$objective)
}
