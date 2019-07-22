# Moment methods for AST & GAT distributions, while gat methods are not yet implemented
# we may also want to keep separate files for both distributions, doesn't seem necessary at the time
# and moment methods for ast & gat class without data, just for exploration uses
# may also want separate functions for mean, variance, sd, skewness & kurtosis
# authorized domain, here or in the plot file

#' @title Asymmetric Student t-distribution
#'
#' @description Probablity density function (pdf), Cumulative distribution function (cdf), quantile function and random generation of AST distributions
#' @param n order of moment to be calculated
#' @param mu location parameter, the mode, not necessarily the mean
#' @param sigma scale parameter, not necessarily the standard deviation, greater than 0
#' @param alpha skewness parameter, ranges from 0 to 1, when < 0.5, skew to the right, when > 0.5, skew to the left
#' @param nu1 degrees of freedom / tail parameter for the left tail, greater than 0
#' @param nu2 degrees of freedom / tail parameter for the right tail, greater than 0
#' @param pars a vector that contains c(mu, sigma, alpha, nu1, nu2), one and only one of {mu, sigma, alpha, nu1, nu2} or pars should be specified
#'
#' @details the moments for AST follow the general rule of student t distribution, mean is only defined for nu > 1, variance is finite when nu > 2, infinite when 1 < nu < 2, otherwise undefined; skewness is defined when nu > 3; kurtosis is finite when nu > 4, infinite when 2 < nu <= 4, otherwise undefined (This need more check, and also the analytical methods should be modified)
#'
#'
#' @aliases moment_ast
#' @aliases mean_ast
#' @aliases var_ast
#' @aliases skew_ast
#' @aliases kurt_ast
#' @aliases moments_ast
#' @name astMoments
#'
#' @examples
#' # The parameter values are specially set for a volatile portfolio.
#' pars <- c(0.12, 0.6, 0.4, 3, 5)
#' moment_ast(1, pars = pars, method = "numerical")
#' mean_ast(pars = pars)
#' var_ast(pars = pars)
#' skew_ast(pars = pars)
#' kurt_ast(pars = pars)
#' moments_ast(pars = pars)

#' @export
moment_ast <- function(n, mu = 0, sigma = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical")) {
  method <- match.arg(method)
  if (!is.null(pars)) {
    if (!missing(mu)) {
      stop("Only one of [mu, sigma, alpha, nu1, nu2] and pars needs to be specified")
    }
    mu <- pars[1]
    sigma <- pars[2]
    alpha <- pars[3]
    nu1 <- pars[4]
    nu2 <- pars[5]
  }
  if (method == "analytical") {
    # return
    mo <- sum( choose(n, 0:n) * sapply(n:0, moment_ss, sigma, alpha, nu1, nu2) * mu^(0:n) )
  }
  if (method == "numerical") {
    integrand <- function(x) {
      x^n * dast(x, mu, sigma, alpha, nu1, nu2)
    }
    # return
    mo <- safeIntegrate(integrand, -Inf, Inf)$value
  }
  mo
}

#' @export
moments_ast <- function(mu = 0, sigma = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical")) {
  if (is.null(pars)) {
    if (!missing(mu)) {
      stop("One and only one of [mu, sigma, alpha, nu1, nu2] and pars needs to be specified")
    }
    pars = c(mu, sigma, alpha, nu1, nu2)
  }
  c(mean = mean_ast(pars = pars, method = method),
    variance = var_ast(pars = pars, method = method),
    skewness = skew_ast(pars = pars, method = method),
    kurtosis = kurt_ast(pars = pars, method = method))
}

#' @export
mean_ast <- function(mu = 0, sigma = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical")) {
  if (!is.null(pars)) {
    if (!missing(mu)) {
      stop("Only one of [mu, sigma, alpha, nu1, nu2] and pars needs to be specified")
    }
    mu <- pars[1]
    sigma <- pars[2]
    alpha <- pars[3]
    nu1 <- pars[4]
    nu2 <- pars[5]
  }
  # return
  moment_ast(1, mu, sigma, alpha, nu1, nu2, method = method)
}

#' @export
var_ast <- function(mu = 0, sigma = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical")) {
  if (!is.null(pars)) {
    if (!missing(mu)) {
      stop("Only one of [mu, sigma, alpha, nu1, nu2] and pars needs to be specified")
    }
    mu <- pars[1]
    sigma <- pars[2]
    alpha <- pars[3]
    nu1 <- pars[4]
    nu2 <- pars[5]
  }
  # return
  moment_central_ast(2, mu, sigma, alpha, nu1, nu2, method = method)
}

#' @export
sd_ast <- function(mu = 0, sigma = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical")) {
  if (!is.null(pars)) {
    if (!missing(mu)) {
      stop("Only one of [mu, sigma, alpha, nu1, nu2] and pars needs to be specified")
    }
    mu <- pars[1]
    sigma <- pars[2]
    alpha <- pars[3]
    nu1 <- pars[4]
    nu2 <- pars[5]
  }
  var <- moment_central_ast(2, mu, sigma, alpha, nu1, nu2, method = method)
  # return
  sqrt(var)
}

#' @export
skew_ast <- function(mu = 0, sigma = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical")) {
  if (!is.null(pars)) {
    if (!missing(mu)) {
      stop("Only one of [mu, sigma, alpha, nu1, nu2] and pars needs to be specified")
    }
    mu <- pars[1]
    sigma <- pars[2]
    alpha <- pars[3]
    nu1 <- pars[4]
    nu2 <- pars[5]
  }
  sd <- sd_ast(mu, sigma, alpha, nu1, nu2, method = method)
  # return
  moment_central_ast(3, mu, sigma, alpha, nu1, nu2, method = method) / sd^3
}

#' @export
kurt_ast <- function(mu = 0, sigma = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical")) {
  if (!is.null(pars)) {
    if (!missing(mu)) {
      stop("Only one of [mu, sigma, alpha, nu1, nu2] and pars needs to be specified")
    }
    mu <- pars[1]
    sigma <- pars[2]
    alpha <- pars[3]
    nu1 <- pars[4]
    nu2 <- pars[5]
  }
  var <- var_ast(mu, sigma, alpha, nu1, nu2, method = method)
  # return
  moment_central_ast(4, mu, sigma, alpha, nu1, nu2, method = method) / var^2
}

moment_central_ast <- function(n, mu, sigma, alpha, nu1, nu2, method = c("analytical", "numerical")) {
  method <- match.arg(method)
  mean <- moment_ast(1, mu, sigma, alpha, nu1, nu2, method = method)
  if (method == "analytical") {
    # return
    mo <- sum( (-1)^(n - n:0) * choose(n, 0:n) * sapply(n:0, moment_ast, mu, sigma, alpha, nu1, nu2) * mean^(0:n) )
  }
  if (method == "numerical") {
    integrand <- function(x) {
      (x - mean)^n * dast(x, mu, sigma, alpha, nu1, nu2)
    }
    # return
    mo <- safeIntegrate(integrand, -Inf, Inf)$value
  }
  mo
}

moment_ss <- function(n, sigma, alpha, nu1, nu2) {
  # moment for sz, s is scale, z is a standardi ast r.v.
  B <- alpha * K(nu1) + (1 - alpha) * K(nu2)
  alpha_star <- alpha * K(nu1)/B
  # return
  alpha * (-2 * alpha_star * sigma * B)^n * moment_abst(nu1, n) +
    (1 - alpha) * (2 * (1 - alpha_star) * sigma * B)^n * moment_abst(nu2, n)
}

moment_abst <- function(nu, n) {
  # absolute moemnt of standard student t
  # -1 < n < nu
  sqrt(nu^n / pi) * gamma( (n+1)/2 ) * gamma( (nu-n)/2 ) / gamma(nu/2)
}
