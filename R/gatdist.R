#' @title Generalized Asymmetric t-distribution
#'
#' @description Probablity density function (pdf), Cumulative distribution function (cdf), quantile function and random generation of GAT distributions
#' @param x,q vector of quantiles
#' @param p vector of probablilities
#' @param n number of observations for random generation
#' @param mu location parameter
#' @param phi scale parameter
#' @param r tail power asymmetry
#' @param c scale asymmetry
#' @param alpha how early tail behavior is apparent
#' @param nu degrees of freedom / tail parameter
#'
#' @aliases dgat
#' @aliases pgat
#' @aliases rgat
#' @name gatDist
#'
#' @examples
#' # The parameter values are specially set for a volatile portfolio.
#' # The sample code are not yet decided because of the incompleteness
#' # of code for GAT distributions.
#' d <- dgat(0, 0.12, 0.6, 0.3, 3, 5)
#' p <- pgat(1.5, 0.12, 0.6, 0.3, 3, 5)
#' q <- qgat(0.8, 0.12, 0.6, 0.3, 3, 5)
#' x <- rgat(1000, 0.12, 0.6, 0.3, 3, 5)



####### all functions needs further checking. Quantile requires newton's method
#' @rdname gatDist
#' @export
dgat <- function(x, mu, phi, alpha, r, c, nu) {
    if (!is.numeric(x))
        stop("x must be numeric")
    g <- (x - mu)/phi + sqrt(1 + ((x - mu)/phi)^2)
    A <- alpha * (1 + r^2)/(r * phi)
    B <- ( (c * g)^(alpha * r) + (c * g)^(-alpha/r) )^(-nu/alpha) / beta(nu/alpha/(1 + r^2), r^2 * nu/alpha/(1 + r^2))
    C <- (1 + ((x - mu)/phi)^2)^(-0.5)
    d <- A * B * C
    d
}

#' @rdname gatDist
#' @export
pgat <- function(x, mu, phi, alpha, r, c, nu) {
    if (!is.numeric(x))
        stop("x must be numeric")
    q <- 1/(1 + c^(-alpha * (1 + r^2)/r) * (((x - mu)/phi) + sqrt(1 + (x - mu)^2/phi^2))^(-alpha * (1 + r^2)/r))
    p <- pbeta(q, nu/alpha/(1 + r^2), r^2 * nu/alpha/(1 + r^2))
    p
}

#' @rdname gatDist
#' @export
qgat <- function(p, mu, phi, alpha, r, c, nu) {
    # tbd, use newton-raphson
  pfunc <- function(q) {
    pgat(q, mu, phi, alpha, r, c, nu) - p
  }
  ps = numeric(21)
  for(i in 1:21) {
    ps[i] <- pfunc(i - 11)
  }
  x0 <- min(abs(ps))
  newton_raphson(pfunc, x0 = x0)
}

#' @rdname gatDist
#' @export
rgat <- function(n, mu, phi, alpha, r, c, nu) {
    if (n < 0)
        stop("x must be non-negative")
    a <- nu/alpha/(1 + r^2)
    b <- nu * r^2/alpha/(1 + r^2)
    q <- rbeta(n, a, b)
    delta <- r/alpha/(1 + r^2)
    x <- mu + 0.5 * phi * ((q/(1 - q))^delta/c - c * (q/(1 - q))^(-delta))
    x
}
