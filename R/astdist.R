#' @title Asymmetric Student t-distribution
#'
#' @description Probablity density function (pdf), Cumulative distribution function (cdf), quantile function and random generation of AST distributions
#' @param x,q vector of quantiles
#' @param p vector of probablilities
#' @param n number of observations for random generation
#' @param mu location parameter, the mode, not necessarily the mean
#' @param sigma scale parameter, not necessarily the standard deviation, greater than 0
#' @param alpha skewness parameter, ranges from 0 to 1, when < 0.5, skew to the right, when > 0.5, skew to the left
#' @param nu1 degrees of freedom / tail parameter for the left tail, greater than 0
#' @param nu2 degrees of freedom / tail parameter for the right tail, greater than 0
#' @param pars a vector that contains c(mu, sigma, alpha, nu1, nu2), one and only one of {mu, sigma, alpha, nu1, nu2} or pars should be specified
#'
#' @aliases dast
#' @aliases past
#' @aliases qast
#' @aliases rast
#' @name astDist
#'
#' @examples
#' # The parameter values are specially set for a volatile portfolio.
#' # density at the mu is always 1 / sigma
#' d <- dast(0.12, 0.12, 0.6, 0.6, 3, 5)
#' # cumulative distribution at mu is alpha
#' p <- past(0.12, 0.12, 0.6, 0.6, 3, 5)
#' # quantile at alpha is mu
#' q <- qast(0.4, 0.12, 0.6, 0.6, 3, 5)
#' data <- rast(1000, 0.12, 0.6, 0.6, 3, 5)
#' hist(data, breaks = 50, probability = TRUE)
#'
#' # now trying to use the pars argument
#' pars <- c(0.12, 0.6, 0.6, 3, 5)
#' x <- seq(-3, 3, 0.01)
#' y <- dast(x, pars = pars)
#' lines(x, y, col = 4)
#'
#' pars1 <- c(0, 2, 0.3, 5, 5)
#' pars2 <- c(0, 2, 0.5, 5, 5)
#' pars3 <- c(0, 2, 0.7, 5, 5)
#' y1 <- dast(x1, pars = pars1)
#' y2 <- dast(x1, pars = pars2)
#' y3 <- dast(x1, pars = pars3)
#' plot(x, y1, type = "l", main = expression(alpha), xlab = "x", ylab = "density")
#' lines(x, y2, lty = 2)
#' lines(x, y3, lty = 3)
#' abline(v = 0, col = 4, lty = 2)
#' legend(x = "topleft", legend = c("alpha = 0.3", "alpha = 0.5", "alpha = 0.7"), lty = 1:3)

#' @rdname astDist
#' @export
dast <- function(x, mu = 0, sigma = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL) {
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
    if (!is.numeric(x))
      stop("x must be numeric")
    x1 <- x[x <= mu]
    x2 <- x[x > mu]
    d <- numeric(length(x))
    # refer to helper_functions.R for K(.)
    d[x <= mu] <- (1 + ((x1 - mu)/(2 * alpha * sigma * K(nu1)))^2/nu1)^(-0.5 * (nu1 + 1))/sigma
    d[x > mu] <- (1 + ((x2 - mu)/(2 * (1 - alpha) * sigma * K(nu2)))^2/nu2)^(-0.5 * (nu2 + 1))/sigma
    d
}

#' @rdname astDist
#' @export
past <- function(q, mu = 0, sigma = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL) {
    if (!is.numeric(q))
        stop("q must be numeric")
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
    B <- alpha * K(nu1) + (1 - alpha) * K(nu2)
    q <- (q - mu) / (sigma * B)
    # refer to helper_functions.R for K(.)
    alpha_star <- alpha * K(nu1)/ B
    # return value
    2 * alpha * pt(pmin(q, 0)/(2 * alpha_star), nu1) + 2 * (1 - alpha) * (pt(pmax(q, 0)/(2 * (1 - alpha_star)), nu2) - 0.5)
}

#' @rdname astDist
#' @export
qast <- function(p, mu = 0, sigma = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL) {
    if (!is.numeric(p))
        stop("p must be numeric")
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
    # refer to helper_functions.R for K(.)
    B <- alpha * K(nu1) + (1 - alpha) * K(nu2)
    alpha_star <- alpha * K(nu1)/ B
    # return value
    mu + sigma * B * (2 * alpha_star * ( qt(pmin(p, alpha)/(2 * alpha), nu1)) + 2 * (1 - alpha_star) * ( qt((pmax(p, alpha) + 1 - 2 * alpha)/(2 * (1 - alpha)), nu2)))
}

#' @rdname astDist
#' @export
rast <- function(n, mu = 0, sigma = 1, alpha = 0.5, nu1 = Inf, nu2 = Inf, pars = NULL) {
    if (n < 0)
        stop("n must be non-negative")
    # refer to helper_functions.R for K(.)
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
    alpha_star <- alpha * K(nu1)/(alpha * K(nu1) + (1 - alpha) * K(nu2))
    u <- runif(n)
    t1 <- rt(n, nu1)
    t2 <- rt(n, nu2)
    x <- alpha_star * abs(t1) * (sign(u - alpha) - 1) + (1 - alpha_star) * abs(t2) * (sign(u - alpha) + 1)
    # return value
    # important note here!  sigma has been transformed
    mu + x * sigma * (alpha * K(nu1) + (1 - alpha) * K(nu2))
}
