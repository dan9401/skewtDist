# log-likelihood function of the AST distributions pars: parameter values y: data which you fit the distribution on
llast <- function(pars, arglist) {
  y <- arglist$data
  ipars <- arglist$ipars

  fixed_pars <- ipars$fixed_pars
  names(fixed_pars) <- rownames(ipars)
  # rename the vector, because nloptr passed vector will remove names
  # if we could only add an assert here, would force the name check, not 100% sure for now
  names(pars) <- rownames(ipars)[which(is.na(fixed_pars))]
  all_pars <- c(pars, fixed_pars[which(!is.na(fixed_pars))])

  mu <- all_pars["mu"]
  sigma <- all_pars["sigma"]
  alpha <- all_pars["alpha"]
  nu1 <- all_pars["nu1"]
  nu2 <- all_pars["nu2"]
  T_ <- length(y)
  y1 <- y[y <= mu]
  y2 <- y[y > mu]

  logl <- -T_ * log(sigma) - 0.5 * (nu1 + 1) * sum(log(1 + ((y1 - mu)/(2 * alpha * sigma * K(nu1)))^2/nu1)) -
    0.5 * (nu2 + 1) * sum(log(1 + ((y2 - mu)/(2 * (1 - alpha) * sigma * K(nu2)))^2/nu2))
  -logl
}

llast_grad <- function(pars, arglist) {
  y <- arglist$data
  ipars <- arglist$ipars

  fixed_pars <- ipars$fixed_pars
  names(fixed_pars) <- rownames(ipars)
  # rename the vector, because nloptr passed vector will remove names
  # if we could only add an assert here, would force the name check, not 100% sure for now
  names(pars) <- rownames(ipars)[which(is.na(fixed_pars))]
  all_pars <- c(pars, fixed_pars[which(!is.na(fixed_pars))])

  mu <- all_pars["mu"]
  sigma <- all_pars["sigma"]
  alpha <- all_pars["alpha"]
  nu1 <- all_pars["nu1"]
  nu2 <- all_pars["nu2"]
  T_ <- length(y)
  y1 <- y[y <= mu]
  y2 <- y[y > mu]

  g_mu <- sum((nu1 + 1) / L(pars, y1) / nu1 * (y1 - mu) / (2 * alpha * sigma * K(nu1))^2) +
    sum((nu2 + 1) / R(pars, y2) / nu2 * (y2 - mu) / (2 * (1 - alpha) * sigma * K(nu2))^2)
  g_sigma<- -T_ / sigma + (nu1 + 1) / sigma * sum(1 - 1 / L(pars, y1)) + (nu2 + 1) / sigma * sum(1 - 1 / R(pars, y2))
  g_alpha <- (nu1 + 1) / alpha * sum(1 - 1 / L(pars, y1)) - (nu2 + 1) / (1 - alpha) * sum(1 - 1 / R(pars, y2))
  g_nu1 <- sum(- log(L(pars, y1)) / 2 + (nu1 + 1) / 2 * D(nu1) * (L(pars, y1) - 1) / L(pars, y1))
  g_nu2 <- sum(- log(R(pars, y2)) / 2 + (nu2 + 1) / 2 * D(nu2) * (R(pars, y2) - 1) / R(pars, y2))
  gradient <- -c(mu = g_mu, sigma = g_sigma, alpha = g_alpha, nu1 = g_nu1, nu2 = g_nu2)

  return(gradient[names(pars)])
}

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
  mu <- pars["mu"]
  alpha <- pars["alpha"]
  nu1 <- pars["nu1"]
  L <- 1 + 1 / nu1 * ( (y - mu) / (2 * alpha * K(nu1)) )^2
  L
}

R <- function(pars, y) {
  mu <- pars["mu"]
  alpha <- pars["alpha"]
  nu2 <- pars["nu2"]
  R <- 1 + 1 / nu2 * ( (y - mu) / (2 * (1 - alpha) * K(nu2)) )^2
  R
}

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
