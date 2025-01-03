#' @title Asymmetric Student-t Distribution
#'
#' @name AST
#' @aliases ast
#' @aliases dast
#' @aliases past
#' @aliases qast
#' @aliases rast
#'
#' @description Probablity density function(PDF), Cumulative distribution function(CDF), Quantile function and Random generation of the AST distribution
#'
#' @param x,q vector of quantiles
#' @param p vector of probablilities
#' @param n number of observations for random generation
#' @param mu location parameter
#' @param s scale parameter, \eqn{s > 0}
#' @param alpha skewness parameter, \eqn{0 < alpha < 1}
#' @param nu1 degrees of freedom / tail parameter for the left tail, \eqn{ nu1 > 0}
#' @param nu2 degrees of freedom / tail parameter for the right tail, \eqn{ nu2 > 0}
#' @param pars a vector that contains mu, s, alpha, nu1, nu2, if pars is specified, mu, s, alpha, nu1, nu2 should not be specified
#'
#' @return
#' \code{dast} gives the density, \code{past} gives the distribution function, \code{qast} gives the quantile function, and \code{rast} generates random samples for AST distribution.
#'
#' @details
#' The 'asymmetric' in AST distribution, not only suggests skewness in the distribution,
#' but also the asymmetry in the two tail powers of the distribution.
#' \itemize{
#'     \item Location parameter \code{mu} is the mode, but not necessarily the mean of the distribution.
#'     \item Scale parameter \code{s} is not necessarily the standard deviation.
#'     \item The distribution skews to the right when the skewness parameter \eqn{alpha < 0.5},
#'     skews to the left when \eqn{alpha > 0.5}.
#'     \item The location paramter \code{mu} always locates at the \eqn{\alpha}-th percentile of the distribution.
#'     The two degrees of freedom / tail parameters each controls one tail of the distribution,
#'     separated at the location paramter \code{mu}.
#'     The left tail parameter \code{nu1} only affects the left half(0th to \eqn{\alpha}-th percentile) of the distribution,
#'     while the right tail paramter \code{nu2} only affects the right half(\eqn{\alpha}-th to 100-th percentile) of the distribution.
#' }
#'
#' @references
#' Zhu, D., & Galbraith, J. W. (2010). A generalized asymmetric Student-t distribution with application to financial econometrics. Journal of Econometrics, 157(2), 297-305.\url{https://www.sciencedirect.com/science/article/pii/S0304407610000266}
#' \url{https://econpapers.repec.org/paper/circirwor/2009s-13.htm}
#'
#' @examples
#' # The parameter values are specially set for a volatile portfolio.
#' # density at the mu is always 1 / s
#' d <- dast(0.12, 0.12, 0.6, 0.6, 3, 5)
#' # cumulative distribution at mu is alpha
#' p <- past(0.12, 0.12, 0.6, 0.6, 3, 5)
#' # quantile at alpha is mu
#' q <- qast(0.4, 0.12, 0.6, 0.6, 3, 5)
#' data <- rast(1000, 0.12, 0.6, 0.6, 3, 5)
#' hist(data, breaks = 50, probability = TRUE)
#'
#' # using the 'pars' argument
#' pars <- c(0.12, 0.6, 0.6, 3, 5)
#' x <- seq(-3, 3, 0.01)
#' y <- dast(x, pars = pars)
#' lines(x, y, col = 4)
#'
#' pars1 <- c(0, 2, 0.3, 5, 5)
#' pars2 <- c(0, 2, 0.5, 5, 5)
#' pars3 <- c(0, 2, 0.7, 5, 5)
#' y1 <- dast(x, pars = pars1)
#' y2 <- dast(x, pars = pars2)
#' y3 <- dast(x, pars = pars3)
#' plot(x, y1, type = "l", main = expression(alpha), xlab = "x", ylab = "density")
#' lines(x, y2, lty = 2)
#' lines(x, y3, lty = 3)
#' abline(v = 0, col = 4, lty = 2)
#' legend(x = "topleft", legend = c("alpha = 0.3", "alpha = 0.5", "alpha = 0.7"), lty = 1:3)
#' 
#' @importFrom stats pt qt rt runif

#' @rdname AST
#' @export
dast <- function(x, mu = 0, s = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL) {
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
    if (!is.numeric(x))
      stop("x must be numeric")
    x1 <- x[x <= mu]
    x2 <- x[x > mu]
    d <- numeric(length(x))
    # refer to helper_functions.R for K(.)
    d[x <= mu] <- (1 + ((x1 - mu)/(2 * alpha * s * K(nu1)))^2/nu1)^(-0.5 * (nu1 + 1))/s
    d[x > mu] <- (1 + ((x2 - mu)/(2 * (1 - alpha) * s * K(nu2)))^2/nu2)^(-0.5 * (nu2 + 1))/s
    d
}

#' @rdname AST
#' @export
past <- function(q, mu = 0, s = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL) {
    if (!is.numeric(q))
        stop("q must be numeric")
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
    B <- alpha * K(nu1) + (1 - alpha) * K(nu2)
    q <- (q - mu) / (s * B)
    # refer to helper_functions.R for K(.)
    alpha_star <- alpha * K(nu1)/ B
    # return value
    2 * alpha * pt(pmin(q, 0)/(2 * alpha_star), nu1) + 2 * (1 - alpha) * (pt(pmax(q, 0)/(2 * (1 - alpha_star)), nu2) - 0.5)
}

#' @rdname AST
#' @export
qast <- function(p, mu = 0, s = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL) {
    if (!is.numeric(p))
        stop("p must be numeric")
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
    # refer to helper_functions.R for K(.)
    B <- alpha * K(nu1) + (1 - alpha) * K(nu2)
    alpha_star <- alpha * K(nu1)/ B
    # return value
    mu + s * B * (2 * alpha_star * ( qt(pmin(p, alpha)/(2 * alpha), nu1)) + 2 * (1 - alpha_star) * ( qt((pmax(p, alpha) + 1 - 2 * alpha)/(2 * (1 - alpha)), nu2)))
}

#' @rdname AST
#' @export
rast <- function(n, mu = 0, s = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL) {
    if (n < 0)
        stop("n must be non-negative")
    # refer to helper_functions.R for K(.)
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
    alpha_star <- alpha * K(nu1)/(alpha * K(nu1) + (1 - alpha) * K(nu2))
    u <- runif(n)
    t1 <- rt(n, nu1)
    t2 <- rt(n, nu2)
    x <- alpha_star * abs(t1) * (sign(u - alpha) - 1) + (1 - alpha_star) * abs(t2) * (sign(u - alpha) + 1)
    # return value
    # important note here!  s has been transformed
    mu + x * s * (alpha * K(nu1) + (1 - alpha) * K(nu2))
}
