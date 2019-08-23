#' @title Generalized Asymmetric t-distribution
#'
#' @name GAT
#' @aliases gat
#' @aliases dgat
#' @aliases pgat
#' @aliases qgat
#' @aliases rgat
#'
#' @description Probablity density function(PDF), Cumulative distribution function(CDF), Quantile function and Random generation of the gat distribution
#'
#' @param x,q vector of quantiles
#' @param p vector of probablilities
#' @param n number of observations for random generation
#' @param mu location parameter
#' @param mu location parameter
#' @param phi scale parameter
#' @param r tail power asymmetry
#' @param c scale asymmetry
#' @param alpha how early tail behavior is apparent
#' @param nu degrees of freedom / tail parameter
#' @param sigma scale parameter, \eqn{sigma > 0}
#' @param alpha skewness parameter, \eqn{0 < alpha < 1}
#' @param nu1 degrees of freedom / tail parameter for the left tail, \eqn{ nu1 > 0}
#' @param nu2 degrees of freedom / tail parameter for the right tail, \eqn{ nu2 > 0}
#' @param pars a vector that contains mu, sigma, alpha, nu1, nu2, if pars is specified, mu, sigma, alpha, nu1, nu2 should not be specified
#'
#' @return
#' \code{dgat} gives the density, \code{pgat} gives the distribution function, \code{qgat} gives the quantile function, and \code{rgat} generates random samples for gat distribution.
#'
#' @details
#' The 'asymmetric' in gat distribution, not only suggests skewness in the distribution,
#' but also the asymmetry in the two tail powers of the distribution.
#' \itemize{
#'     \item Location parameter \code{mu} is the mode, but not necessarily the mean of the distribution.
#'     \item Scale parameter \code{sigma} is not necessarily the standard deviation.
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
#' # density at the mu is always 1 / sigma
#' d <- dgat(0.12, 0.12, 0.6, 0.6, 3, 5)
#' # cumulative distribution at mu is alpha
#' p <- pgat(0.12, 0.12, 0.6, 0.6, 3, 5)
#' # quantile at alpha is mu
#' q <- qgat(0.4, 0.12, 0.6, 0.6, 3, 5)
#' data <- rgat(1000, 0.12, 0.6, 0.6, 3, 5)
#' hist(data, breaks = 50, probability = TRUE)
#'
#' # using the 'pars' argument
#' pars <- c(0.12, 0.6, 0.6, 3, 5)
#' x <- seq(-3, 3, 0.01)
#' y <- dgat(x, pars = pars)
#' lines(x, y, col = 4)
#'
#' pars1 <- c(0, 2, 0.3, 5, 5)
#' pars2 <- c(0, 2, 0.5, 5, 5)
#' pars3 <- c(0, 2, 0.7, 5, 5)
#' y1 <- dgat(x, pars = pars1)
#' y2 <- dgat(x, pars = pars2)
#' y3 <- dgat(x, pars = pars3)
#' plot(x, y1, type = "l", main = expression(alpha), xlab = "x", ylab = "density")
#' lines(x, y2, lty = 2)
#' lines(x, y3, lty = 3)
#' abline(v = 0, col = 4, lty = 2)
#' legend(x = "topleft", legend = c("alpha = 0.3", "alpha = 0.5", "alpha = 0.7"), lty = 1:3)

#' @rdname GAT
#' @export
dgat <- function(x, mu, phi, alpha, r, c, nu, pars = NULL) {
    if (!is.numeric(x))
        stop("x must be numeric")
    if (!is.null(pars)) {
      if (!missing(mu)) {
        stop("Only one of [mu, sigma, alpha, nu1, nu2] and pars needs to be specified")
      }
      mu <- pars[1]
      phi <- pars[2]
      alpha <- pars[3]
      r <- pars[4]
      c <- pars[5]
      nu <- pars[6]
    }
    g <- (x - mu)/phi + sqrt(1 + ((x - mu)/phi)^2)
    A <- alpha * (1 + r^2)/(r * phi)
    B <- ( (c * g)^(alpha * r) + (c * g)^(-alpha/r) )^(-nu/alpha) / beta(nu/alpha/(1 + r^2), r^2 * nu/alpha/(1 + r^2))
    C <- (1 + ((x - mu)/phi)^2)^(-0.5)
    d <- A * B * C
    d
}

#' @rdname GAT
#' @export
pgat <- function(x, mu, phi, alpha, r, c, nu, pars = NULL) {
    if (!is.numeric(x))
        stop("x must be numeric")
    if (!is.null(pars)) {
      if (!missing(mu)) {
        stop("Only one of [mu, sigma, alpha, nu1, nu2] and pars needs to be specified")
      }
      mu <- pars[1]
      phi <- pars[2]
      alpha <- pars[3]
      r <- pars[4]
      c <- pars[5]
      nu <- pars[6]
    }
    if (x == -Inf) {
      return(0)
    }
    if (x == Inf) {
      return(1)
    }
    q <- 1/(1 + c^(-alpha * (1 + r^2)/r) * (((x - mu)/phi) + sqrt(1 + (x - mu)^2/phi^2))^(-alpha * (1 + r^2)/r))
    p <- pbeta(q, nu/alpha/(1 + r^2), r^2 * nu/alpha/(1 + r^2))
    p
}

#' @rdname GAT
#' @export
qgat <- function(p, mu, phi, alpha, r, c, nu, pars = NULL) {
    # tbd, use newton-raphson
    if (!is.null(pars)) {
      if (!missing(mu)) {
        stop("Only one of [mu, sigma, alpha, nu1, nu2] and pars needs to be specified")
      }
      mu <- pars[1]
      phi <- pars[2]
      alpha <- pars[3]
      r <- pars[4]
      c <- pars[5]
      nu <- pars[6]
    }
    n <- length(p)
    ps = matrix(0, nrow = 21, ncol = n)

    f <- function(y) {
      uniroot(function(x) {
        pgat(x, mu, phi, alpha, r, c, nu) - y
        },
        interval = c(-1e5, 1e5), tol = 1e-8)$root
    }

    q <- sapply(p, f)
    q
}

#' @rdname GAT
#' @export
rgat <- function(n, mu, phi, alpha, r, c, nu, pars = NULL) {
    if (!is.null(pars)) {
      if (!missing(mu)) {
        stop("Only one of [mu, sigma, alpha, nu1, nu2] and pars needs to be specified")
      }
      mu <- pars[1]
      phi <- pars[2]
      alpha <- pars[3]
      r <- pars[4]
      c <- pars[5]
      nu <- pars[6]
    }
    if (n < 0)
        stop("x must be non-negative")
    a <- nu/alpha/(1 + r^2)
    b <- nu * r^2/alpha/(1 + r^2)
    q <- rbeta(n, a, b)
    delta <- r/alpha/(1 + r^2)
    x <- mu + 0.5 * phi * ((q/(1 - q))^delta/c - c * (q/(1 - q))^(-delta))
    x
}
