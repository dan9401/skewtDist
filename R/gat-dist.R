#' @title Generalized Asymmetric t-distribution
#'
#' @name GAT
#' @aliases gat
#' @aliases dgat
#' @aliases pgat
#' @aliases qgat
#' @aliases rgat
#'
#' @description Probablity density function(PDF), Cumulative distribution function(CDF), Quantile function and Random generation of the GAT distribution
#'
#' @param x,q vector of quantiles
#' @param p vector of probablilities
#' @param n number of observations for random generation
#' @param mu location parameter
#' @param scale scale parameter, \eqn{scale > 0}
#' @param alpha how early tail behavior is apparent, \eqn{alpha > 0}
#' @param r tail power asymmetry, \eqn{r > 0}
#' @param c scale asymmetry, \eqn{c > 0}
#' @param nu degrees of freedom / tail parameter, \eqn{nu > 0}
#' @param pars a vector that contains mu, scale, alpha, r, c, nu, if pars is specified, mu, scale, alpha, r, c, nu should not be specified
#'
#' @return
#' \code{dgat} gives the density, \code{pgat} gives the distribution function, \code{qgat} gives the quantile function, and \code{rgat} generates random samples for GATdistribution.
#'
#' @examples
#' dgat(0, 0, 1, 1.2, 1.2, 2, 5)
#' # using the 'pars' argument
#' pars <- c(0, 1, 1.2, 1.2, 2, 5)
#' x <- seq(-3, 3, 0.01)
#' y <- dgat(x, pars = pars)
#' lines(x, y, col = 4)
#' 
#' @importFrom stats pbeta rbeta uniroot

#' @rdname GAT
#' @export
dgat <- function(x, mu, scale, alpha, r, c, nu, pars = NULL) {
    if (!is.numeric(x))
        stop("x must be numeric")
    if (!is.null(pars)) {
      if (!missing(mu)) {
        stop("Only one of [mu, scale, alpha, r, c, nu] and pars needs to be specified")
      }
      mu <- pars[1]
      scale <- pars[2]
      alpha <- pars[3]
      r <- pars[4]
      c <- pars[5]
      nu <- pars[6]
    }
    g <- (x - mu)/scale + sqrt(1 + ((x - mu)/scale)^2)
    A <- alpha * (1 + r^2)/(r * scale)
    B <- ( (c * g)^(alpha * r) + (c * g)^(-alpha/r) )^(-nu/alpha) / beta(nu/alpha/(1 + r^2), r^2 * nu/alpha/(1 + r^2))
    C <- (1 + ((x - mu)/scale)^2)^(-0.5)
    d <- A * B * C
    d
}

#' @rdname GAT
#' @export
pgat <- function(x, mu, scale, alpha, r, c, nu, pars = NULL) {
    if (!is.numeric(x))
        stop("x must be numeric")
    if (!is.null(pars)) {
      if (!missing(mu)) {
        stop("Only one of [mu, scale, alpha, r, c, nu] and pars needs to be specified")
      }
      mu <- pars[1]
      scale <- pars[2]
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
    q <- 1/(1 + c^(-alpha * (1 + r^2)/r) * (((x - mu)/scale) + sqrt(1 + (x - mu)^2/scale^2))^(-alpha * (1 + r^2)/r))
    p <- pbeta(q, nu/alpha/(1 + r^2), r^2 * nu/alpha/(1 + r^2))
    p
}

#' @rdname GAT
#' @export
qgat <- function(p, mu, scale, alpha, r, c, nu, pars = NULL) {
    # tbd, use newton-raphson
    if (!is.null(pars)) {
      if (!missing(mu)) {
        stop("Only one of [mu, scale, alpha, r, c, nu] and pars needs to be specified")
      }
      mu <- pars[1]
      scale <- pars[2]
      alpha <- pars[3]
      r <- pars[4]
      c <- pars[5]
      nu <- pars[6]
    }
    n <- length(p)
    ps = matrix(0, nrow = 21, ncol = n)

    f <- function(y) {
      uniroot(function(x) {
        pgat(x, mu, scale, alpha, r, c, nu) - y
        },
        interval = c(-1e5, 1e5), tol = 1e-8)$root
    }

    q <- sapply(p, f)
    q
}

#' @rdname GAT
#' @export
rgat <- function(n, mu, scale, alpha, r, c, nu, pars = NULL) {
    if (!is.null(pars)) {
      if (!missing(mu)) {
        stop("Only one of [mu, scale, alpha, r, c, nu] and pars needs to be specified")
      }
      mu <- pars[1]
      scale <- pars[2]
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
    x <- mu + 0.5 * scale * ((q/(1 - q))^delta/c - c * (q/(1 - q))^(-delta))
    x
}
