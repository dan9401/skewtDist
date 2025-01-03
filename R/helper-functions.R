#' @importFrom stats dnorm

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
  Dp <- 1/2 * (trigamma((nu+1)/2) - trigamma(nu/2))
  Dp
}

L <- function(pars, y) {
  mu <- pars[1]
  sigma <- pars[2]
  alpha <- pars[3]
  nu1 <- pars[4]
  L <- 1 + 1 / nu1 * ( (y - mu) / (2 * alpha * sigma * K(nu1)) )^2
  L
}

R <- function(pars, y) {
  mu <- pars[1]
  sigma <- pars[2]
  alpha <- pars[3]
  nu2 <- ifelse(length(pars) == 4, pars[4], pars[5])
  R <- 1 + 1 / nu2 * ( (y - mu) / (2 * (1 - alpha) * sigma * K(nu2)) )^2
  R
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

# putting it here temporarily
#' @export
moments <- function(x, ...) {
  UseMethod("moments", x)
}
