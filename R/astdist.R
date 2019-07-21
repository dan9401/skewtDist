#' @title Asymmetric Student t-distribution
#'
#' @description Probablity density function (pdf), Cumulative distribution function (cdf), quantile function and random generation of AST distributions
#' @param x,q vector of quantiles
#' @param p vector of probablilities
#' @param n number of observations for random generation
#' @param mu location parameter
#' @param sigma scale parameter
#' @param alpha skewness parameter
#' @param nu1 degrees of freedom / tail parameter for the left tail
#' @param nu2 degrees of freedom / tail parameter for the right tail
#'
#' @aliases dast
#' @aliases past
#' @aliases qast
#' @aliases rast
#' @name astDist
#'
#' @examples
#' # The parameter values are specially set for a volatile portfolio.
#' d <- dast(0, 0.12, 0.6, 0.3, 3, 5)
#' p <- past(1.5, 0.12, 0.6, 0.3, 3, 5)
#' q <- qast(0.8, 0.12, 0.6, 0.3, 3, 5)
#' x <- rast(1000, 0.12, 0.6, 0.3, 3, 5)

#' @rdname astDist
#' @export
dast.numeric <- function(x, mu, sigma, alpha, nu1, nu2) {
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
past.numeric <- function(q, mu, sigma, alpha, nu1, nu2) {
    if (!is.numeric(q))
        stop("q must be numeric")
    q <- (q - mu)/sigma
    # refer to helper_functions.R for K(.)
    alpha_star <- alpha * K(nu1)/(alpha * K(nu1) + (1 - alpha) * K(nu2))
    # return value
    2 * alpha * pt(pmin(q, 0)/(2 * alpha_star), nu1) + 2 * (1 - alpha) * (pt(pmax(q, 0)/(2 * (1 - alpha_star)), nu2) - 0.5)
}

#' @rdname astDist
#' @export
qast.numeric <- function(p, mu, sigma, alpha, nu1, nu2) {
    if (!is.numeric(p))
        stop("p must be numeric")
    # refer to helper_functions.R for K(.)
    alpha_star <- alpha * K(nu1)/(alpha * K(nu1) + (1 - alpha) * K(nu2))
    # return value
    2 * alpha_star * (mu + sigma * qt(pmin(p, alpha)/(2 * alpha), nu1)) + 2 * (1 - alpha_star) * (mu + sigma * qt((pmax(p, alpha) + 1 - 2 * alpha)/(2 * (1 - alpha)), nu2))
}

#' @rdname astDist
#' @export
rast.numeric <- function(n, mu, sigma, alpha, nu1, nu2) {
    if (n < 0)
        stop("n must be non-negative")
    # refer to helper_functions.R for K(.)
    alpha_star <- alpha * K(nu1)/(alpha * K(nu1) + (1 - alpha) * K(nu2))
    u <- runif(n)
    t1 <- rt(n, nu1)
    t2 <- rt(n, nu2)
    x <- alpha_star * abs(t1) * (sign(u - alpha) - 1) + (1 - alpha_star) * abs(t2) * (sign(u - alpha) + 1)
    # return value
    # important note here!  sigma has been transformed
    mu + x * sigma * (alpha * K(nu1) + (1 - alpha) * K(nu2))
}
