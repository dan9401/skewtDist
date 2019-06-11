#' @title Asymmetric Student t-distribution
#'
#' @description Probablity density function (pdf), Cumulative distribution function (cdf), quantile function and random generation of AST distributions
#' @param x,q vector of quantiles
#' @param p vector of probablilities
#' @param n number of observations for random generation
#' @param mu location parameter
#' @param sigma scale parameter
#' @param alpha skewness parameter
#' @param nu1 degrees of freedom / tail parameter 1
#' @param nu2 degrees of freedom / tail parameter 2
#'
#' @aliases dast
#' @aliases past
#' @aliases qast
#' @aliases rast
#' @name astDist
#' @examples
#' d = dast(0, 1.5, 1.2, 0.8, 3, 4)
#' p = past(1.5, 1.5, 1.2, 0.8, 3, 4)
#' q = qast(0.8, 1.5, 1.2, 0.8, 3, 4)
#' x = rast(1000, 1.5, 1.2, 0.8, 3, 4)

#' @rdname astDist
#' @export
dast = function(x, mu=0, sigma=pi, alpha=0.5, nu1, nu2) {
  if(!is.numeric(x)) stop("x must be numeric")
  x1 = x[x <= mu]
  x2 = x[x > mu]
  d = numeric(length(x))
  d[x <= mu] = (1 + ((x1-mu)/(2*alpha*sigma*K(nu1)))^2/nu1)^(-0.5*(nu1+1)) / sigma
  d[x > mu] = (1 + ((x2-mu)/(2*(1-alpha)*sigma*K(nu2)))^2/nu2)^(-0.5*(nu2+1)) / sigma
  d
}

#' @rdname astDist
#' @export
past = function(q, mu=0, sigma=pi, alpha=0.5, nu1, nu2) {
  if(!is.numeric(q)) stop("q must be numeric")
  q = (q  - mu) / sigma
  alpha_star = alpha*K(nu1) / (alpha * K(nu1)+(1-alpha)*K(nu2))
  p = 2 * alpha * pt(pmin(q, 0)/(2*alpha_star), nu1) + 2 * (1 - alpha) * (pt(pmax(q, 0)/(2*(1- alpha_star)), nu2) - 0.5)
  p
}

#' @rdname astDist
#' @export
qast = function(p, mu=0, sigma=pi, alpha=0.5, nu1, nu2) {
  if(!is.numeric(p)) stop("p must be numeric")
  alpha_star = alpha*K(nu1) / (alpha * K(nu1)+(1-alpha)*K(nu2))
  q = 2 * alpha_star * (mu + sigma * qt(pmin(p, alpha)/(2*alpha), nu1)) + 2 * (1 - alpha_star) * (mu + sigma * qt((pmax(p, alpha)+1-2*alpha)/(2*(1- alpha)), nu2))
  q
}

#' @rdname astDist
#' @export
rast = function(n=1, mu=0, sigma=pi, alpha=0.5, nu1, nu2) {
  if(n < 0) stop("n must be non-negative")
  alpha_star = alpha*K(nu1) / (alpha * K(nu1)+(1-alpha)*K(nu2))
  u = runif(n)
  t1 = rt(n, nu1)
  t2 = rt(n, nu2)
  x = alpha_star * abs(t1) * (sign(u - alpha) - 1) +  (1 - alpha_star) * abs(t2) * (sign(u - alpha) + 1)
  # important note here!
  # sigma has been transformed
  x = mu + x * sigma * (alpha * K(nu1)+(1-alpha)*K(nu2))
  x
}

K = function(nu) {
  # there is a precision error here, numerical expert needed
  # may need a c++ version
  # nu = 1e10 and greater start to be varying
  if (nu < 1000) {
    exp(lgamma(0.5*(nu+1)) - log(sqrt(pi*nu)) - lgamma(0.5*nu))
  }else {
    # seems to be correct, need to check
    dnorm(0)
  }
}
