#' @title Moment Functions of Asymmetric Student-t distribution
#'
#' @name ast-moment
#' @aliases astMean
#' @aliases astVar
#' @aliases astSkew
#' @aliases astKurt
#' @aliases astMoment
#' @aliases astMoments
#' @aliases astRawMoment
#' @aliases astCentralMoment
#'
#' @description The mean, standard deviation, skewness, kurtosis functions, as well as the raw and central moments of AST distribution
#'
#' @param moment the moment to be calculated, one of 'mean', 'sd', 'skew', 'kurt'
#' @param n order of (raw/central) moment to be calculated
#' @param mu location parameter
#' @param s scale parameter, \eqn{s > 0}
#' @param alpha skewness parameter, \eqn{0 < alpha < 1}
#' @param nu1 degrees of freedom / tail parameter for the left tail, \eqn{ nu1 > 0}
#' @param nu2 degrees of freedom / tail parameter for the right tail, \eqn{ nu2 > 0}
#' @param pars a vector that contains mu, s, alpha, nu1, nu2, if pars is specified, mu, s, alpha, nu1, nu2 should not be specified
#' @param method method used to calculate the moment(s), one of 'analytical' and 'numerical'
#' @param type type of kurtosis calculated, one of 'excess' and 'regular'
#'
#' @details
#' Function \code{astMoment} calculates one of mean, standard deviation, skewness and kurtosis of the distribution,
#' while \code{astMoment} calculates all 4 of them. \cr
#' Function \code{astRawMoment} returns \eqn{E[X^n]},
#' while function \code{astCentralMoment} returns \eqn{E[(X-\mu)^n]}
#'
#' The moments for AST follow the general rule of student t distribution,
#' \itemize{
#'     \item mean is only defined for nu > 1,
#'     \item variance/standard deviation is finite when nu > 2, infinite when 1 < nu < 2, otherwise undefined,
#'     \item skewness is defined when nu > 3,
#'     \item kurtosis is finite when nu > 4, infinite when 2 < nu <= 4, otherwise undefined.
#' }
#'
#' @references
#' Zhu, D., & Galbraith, J. W. (2010). A generalized asymmetric Student-t distribution with application to financial econometrics. Journal of Econometrics, 157(2), 297-305.\url{https://www.sciencedirect.com/science/article/pii/S0304407610000266}
#' \url{https://econpapers.repec.org/paper/circirwor/2009s-13.htm}
#'
#' @examples
#' # The parameter values are specially set for a volatile portfolio.
#' pars <- c(0.12, 0.6, 0.6, 6, 5)
#' astMoment("sd", pars = pars, method = "numerical")
#' astMoments(pars = pars)
#' 
#' @importFrom DistributionUtils safeIntegrate

#' @export
astMean <- function(mu = 0, s = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical")) {
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, s, alpha, nu1, nu2] and pars needs to be specified")
    }
    pars = c(mu, s, alpha, nu1, nu2)
  }
  # return
  astRawMoment(1, pars = pars, method = method)
}

#' @export
astVar <- function(mu = 0, s = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical")) {
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, s, alpha, nu1, nu2] and pars needs to be specified")
    }
    pars = c(mu, s, alpha, nu1, nu2)
  }
  # return
  astCentralMoment(2, pars = pars, method = method)
}

#' @export
astSD <- function(mu = 0, s = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical")) {
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, s, alpha, nu1, nu2] and pars needs to be specified")
    }
    pars = c(mu, s, alpha, nu1, nu2)
  }
  var <- astCentralMoment(2, pars = pars, method = method)
  # return
  sqrt(var)
}

#' @export
astSkew <- function(mu = 0, s = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical")) {
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, s, alpha, nu1, nu2] and pars needs to be specified")
    }
    pars = c(mu, s, alpha, nu1, nu2)
  }
  sd <- astSD(pars = pars, method = method)
  # return
  astCentralMoment(3, pars = pars, method = method) / sd^3
}

#' @export
astKurt <- function(mu = 0, s = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical"), type = c("excess", "regular")) {
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, s, alpha, nu1, nu2] and pars needs to be specified")
    }
    pars = c(mu, s, alpha, nu1, nu2)
  }
  var <- astVar(pars = pars, method = method)
  kurt <- astCentralMoment(4, pars = pars, method = method) / var^2
  type <- match.arg(type)
  # return
  switch(type,
         excess = kurt - 3,
         regular = kurt)
}

#' @rdname ast-moment
#' @export
astMoment <- function(moment = c("mean", "sd", "var", "skew", "kurt"), mu = 0, s = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical"), type = c("excess", "regular")) {
  moment <- match.arg(moment)
  method <- match.arg(method)
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, s, alpha, nu1, nu2] and pars needs to be specified")
    }
    pars = c(mu, s, alpha, nu1, nu2)
  }
  switch(moment,
         mean = astMean(pars = pars, method = method),
         sd = astSD(pars = pars, method = method),
         var = astVar(pars = pars, method = method),
         skew = astSkew(pars = pars, method = method),
         kurt = astKurt(pars = pars, method = method, type = type))
}

#' @rdname ast-moment
#' @export
astMoments <- function(mu = 0, s = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical"), type = c("excess", "regular")) {
  if (is.null(pars)) {
    if (missing(mu)) {
      stop("One and only one of [mu, s, alpha, nu1, nu2] and pars needs to be specified")
    }
    pars = c(mu, s, alpha, nu1, nu2)
  }
  c(mean = astMean(pars = pars, method = method),
    sd = astSD(pars = pars, method = method),
    skew = astSkew(pars = pars, method = method),
    kurt = astKurt(pars = pars, method = method, type = type))
}

#' @rdname ast-moment
#' @export
astRawMoment <- function(n, mu = 0, s = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL, method = c("analytical", "numerical")) {
  method <- match.arg(method)
  if (!is.null(pars)) {
    if (!missing(mu)) {
      stop("Only one of [mu, s, alpha, nu1, nu2] and pars needs to be specified")
    }
    mu <- pars[1]
    s <- pars[2]
    alpha <- pars[3]
    nu1 <- pars[4]
    nu2 <- pars[5]
  }
  if (method == "analytical") {
    # return
    moment <- sum( choose(n, 0:n) * sapply(n:0, scaledStdASTMoment, s, alpha, nu1, nu2) * mu^(0:n) )
  }
  if (method == "numerical") {
    integrand <- function(x) {
      x^n * dast(x, mu, s, alpha, nu1, nu2)
    }
    # return
    moment <- safeIntegrate(integrand, -Inf, Inf)$value
  }
  moment
}

#' @rdname ast-moment
#' @export
astCentralMoment <- function(n, mu, s, alpha, nu1, nu2, pars = NULL, method = c("analytical", "numerical")) {
  method <- match.arg(method)
  if (!is.null(pars)) {
    if (!missing(mu)) {
      stop("Only one of [mu, s, alpha, nu1, nu2] and pars needs to be specified")
    }
    mu <- pars[1]
    s <- pars[2]
    alpha <- pars[3]
    nu1 <- pars[4]
    nu2 <- pars[5]
  }
  mean <- astRawMoment(1, mu, s, alpha, nu1, nu2, method = method)
  if (method == "analytical") {
    # return
    moment <- sum( (-1)^(n - n:0) * choose(n, 0:n) * sapply(n:0, astRawMoment, mu, s, alpha, nu1, nu2) * mean^(0:n) )
  }
  if (method == "numerical") {
    integrand <- function(x) {
      (x - mean)^n * dast(x, mu, s, alpha, nu1, nu2)
    }
    # return
    moment <- safeIntegrate(integrand, -Inf, Inf)$value
  }
  moment
}

scaledStdASTMoment <- function(n, s, alpha, nu1, nu2) {
  # moment for sz, s is scale, z is a standardi ast r.v.
  B <- alpha * K(nu1) + (1 - alpha) * K(nu2)
  alpha_star <- alpha * K(nu1)/B
  # return
  alpha * (-2 * alpha_star * s * B)^n * absTMoments(nu1, n) +
    (1 - alpha) * (2 * (1 - alpha_star) * s * B)^n * absTMoments(nu2, n)
}

absTMoments <- function(nu, n) {
  # absolute moemnt of standard student t
  # -1 < n < nu
  sqrt(nu^n / pi) * gamma( (n+1)/2 ) * gamma( (nu-n)/2 ) / gamma(nu/2)
}
