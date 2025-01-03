#' @title Moment Functions of Asymmetric Student-t distribution
#'
#' @name gat-moment
#' @aliases gatMean
#' @aliases gatVar
#' @aliases gatSD
#' @aliases gatSkew
#' @aliases gatKurt
#' @aliases gatMoment
#' @aliases gatMoments
#' @aliases gatRawMoment
#' @aliases gatCentralMoment
#'
#' @description The mean, standard deviation, skewness, kurtosis functions, as well as the raw and central moments of GAT distribution
#'
#' @param moment the moment to be calculated, one of 'mean', 'sd', 'skew', 'kurt'
#' @param n order of (raw/central) moment to be calculated
#' @param mu location parameter
#' @param phi scale parameter, \eqn{phi > 0}
#' @param alpha skewness parameter, \eqn{0 < alpha < 1}
#' @param r tail power asymmetry parameter \eqn{r > 0}
#' @param c scale asymmetry parameter \eqn{r > 0}
#' @param nu degrees of freedom / tail parameter
#' @param pars a vector that contains mu, phi, alpha, r, c, nu, if pars is specified, mu, phi, alpha, r, c, nu should not be specified
#' @param method method used to calculate the moment(s), one of 'analytical' and 'numerical'
#' @param type type of kurtosis calculated, one of 'excess' and 'regular'
#'
#' @details
#' Function \code{gatMoment} calculates one of mean, standard deviation, skewness and kurtosis of the distribution,
#' while \code{gatMoment} calculates all 4 of them. \cr
#' Function \code{gatRawMoment} returns \eqn{E[X^n]},
#' while function \code{gatCentralMoment} returns \eqn{E[(X-\mu)^n]}
#'
#' The moments for GAT follow the general rule of student t distribution,
#' \itemize{
#'     \item mean is only defined for nu > 1,
#'     \item variance/standard deviation is finite when nu > 2, infinite when 1 < nu < 2, otherwise undefined,
#'     \item skewness is defined when nu > 3,
#'     \item kurtosis is finite when nu > 4, infinite when 2 < nu <= 4, otherwise undefined.
#' }
#'
#' @examples
#' # The parameter values are specially set for a volatile portfolio.
#' pars <- c(0.12, 0.6, 0.6, 6, 5)
#' gatMoment("sd", pars = pars, method = "numerical")
#' gatMoments(pars = pars)

#' @export
gatMean <- function(mu = 0, phi = 1, alpha = 0.5, r = 2, c = 2, nu = Inf, pars = NULL, method = c("analytical", "numerical")) {
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, phi, alpha, r, c, nu] and pars needs to be specified")
    }
    pars = c(mu, phi, alpha, r, c, nu)
  }
  # return
  gatRawMoment(1, pars = pars, method = method)
}

#' @export
gatVar <- function(mu = 0, phi = 1, alpha = 0.5, r = 2, c = 2, nu = Inf, pars = NULL, method = c("analytical", "numerical")) {
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, phi, alpha, r, c, nu] and pars needs to be specified")
    }
    pars = c(mu, phi, alpha, r, c, nu)
  }
  # return
  gatCentralMoment(2, pars = pars, method = method)
}

#' @export
gatSD <- function(mu = 0, phi = 1, alpha = 0.5, r = 2, c = 2, nu = Inf, pars = NULL, method = c("analytical", "numerical")) {
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, phi, alpha, r, c, nu] and pars needs to be specified")
    }
    pars = c(mu, phi, alpha, r, c, nu)
  }
  var <- gatCentralMoment(2, pars = pars, method = method)
  # return
  sqrt(var)
}

#' @export
gatSkew <- function(mu = 0, phi = 1, alpha = 0.5, r = 2, c = 2, nu = Inf, pars = NULL, method = c("analytical", "numerical")) {
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, phi, alpha, r, c, nu] and pars needs to be specified")
    }
    pars = c(mu, phi, alpha, r, c, nu)
  }
  sd <- gatSD(pars = pars, method = method)
  # return
  gatCentralMoment(3, pars = pars, method = method) / sd^3
}

#' @export
gatKurt <- function(mu = 0, phi = 1, alpha = 0.5, r = 2, c = 2, nu = Inf, pars = NULL, method = c("analytical", "numerical"), type = c("excess", "regular")) {
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, phi, alpha, r, c, nu] and pars needs to be specified")
    }
    pars = c(mu, phi, alpha, r, c, nu)
  }
  var <- gatVar(pars = pars, method = method)
  kurt <- gatCentralMoment(4, pars = pars, method = method) / var^2
  type <- match.arg(type)
  # return
  switch(type,
         excess = kurt - 3,
         regular = kurt)
}

#' @export
gatMoment <- function(moment = c("mean", "sd", "var", "skew", "kurt"), mu = 0, phi = 1, alpha = 0.5, r = 2, c = 2, nu = Inf, pars = NULL, method = c("analytical", "numerical"), type = c("excess", "regular")) {
  moment <- match.arg(moment)
  method <- match.arg(method)
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, phi, alpha, r, c, nu] and pars needs to be specified")
    }
    pars = c(mu, phi, alpha, r, c, nu)
  }
  switch(moment,
         mean = gatMean(pars = pars, method = method),
         sd = gatSD(pars = pars, method = method),
         var = gatVar(pars = pars, method = method),
         skew = gatSkew(pars = pars, method = method),
         kurt = gatKurt(pars = pars, method = method, type = type))
}

#' @rdname gat-moment
#' @export
gatMoments <- function(mu = 0, phi = 1, alpha = 0.5, r = 2, c = 2, nu = Inf, pars = NULL, method = c("analytical", "numerical"), type = c("excess", "regular")) {
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, phi, alpha, r, c, nu] and pars needs to be specified")
    }
    pars = c(mu, phi, alpha, r, c, nu)
  }
  c(mean = gatMean(pars = pars, method = method),
    sd = gatSD(pars = pars, method = method),
    skew = gatSkew(pars = pars, method = method),
    kurt = gatKurt(pars = pars, method = method, type = type))
}

#' @rdname gat-moment
#' @export
gatRawMoment <- function(n, mu = 0, phi = 1, alpha = 0.5, r = 2, c = 2, nu = Inf, pars = NULL, method = c("analytical", "numerical")) {
  method <- match.arg(method)
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, phi, alpha, r, c, nu] and pars needs to be specified")
    }
    pars = c(mu, phi, alpha, r, c, nu)
  }
  if (method == "analytical") {
    # return
    moment <- sum( choose(n, 0:n) * sapply(n:0, scaledStdGATMoment, pars = pars) * mu^(0:n) )
  }
  if (method == "numerical") {
    integrand <- function(x) {
      x^n * dgat(x, pars = pars)
    }
    # return
    moment <- safeIntegrate(integrand, -Inf, Inf)$value
  }
  moment
}

#' @rdname gat-moment
#' @export
gatCentralMoment <- function(n, mu = 0, phi = 1, alpha = 0.5, r = 2, c = 2, nu = Inf, pars = NULL, method = c("analytical", "numerical")) {
  method <- match.arg(method)
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, phi, alpha, r, c, nu] and pars needs to be specified")
    }
    pars = c(mu, phi, alpha, r, c, nu)
  }
  mean <- gatRawMoment(1, pars = pars, method = method)
  if (method == "analytical") {
    # return
    moment <- sum( (-1)^(n - n:0) * choose(n, 0:n) * sapply(n:0, gatRawMoment, pars = pars) * mean^(0:n) )
  }
  if (method == "numerical") {
    integrand <- function(x) {
      (x - mean)^n * dgat(x, pars = pars)
    }
    # return
    moment <- safeIntegrate(integrand, -Inf, Inf)$value
  }
  moment
}

scaledStdGATMoment <- function(n, phi, alpha, r, c, nu, pars = NULL) {
  if (!is.null(pars)) {
    if (!missing(phi)) {
      stop("One and only one of [mu, phi, alpha, r, c, nu] and pars needs to be specified")
    }
    phi <- pars[2]
    alpha <- pars[3]
    r <- pars[4]
    c <- pars[5]
    nu <- pars[6]
  }
  A <- nu/( alpha*(1+r^2) )
  B <- A * r^2
  delta <- A / nu * r
  m <- 0:n
  (phi/2)^n / beta(A, B) * sum( (-1)^m * choose(n, m) * c^(n - 2*m) *
                                            beta(A - (n-2*m)*delta, B + (n-2*m)*delta))
}
