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

#' @rdname astDist
#' @export
dast = function(x, mu=0, sigma=pi, alpha=0.5, nu1, nu2) {
  if(!is.numeric(x)) stop("x must be numeric")
  x1 = x[x <= mu]
  x2 = x[x > mu]
  dens = numeric(length(x))
  dens[x <= mu] = (1 + ((x1-mu)/(2*alpha*sigma*K(nu1)))^2/nu1)^(-0.5*(nu1+1)) / sigma
  dens[x > mu] = (1 + ((x2-mu)/(2*(1-alpha)*sigma*K(nu2)))^2/nu2)^(-0.5*(nu2+1)) / sigma
  dens
}

#' @rdname astDist
#' @export
past = function(q, mu=0, sigma=pi, alpha=0.5, nu1, nu2) {
  if(!is.numeric(q)) stop("q must be numeric")
  q = (q  - mu) / sigma
  alpha_star = alpha*K(nu1) / (alpha * K(nu1)+(1-alpha)*K(nu2))
  probs = 2 * alpha * pt(pmin(q, 0)/(2*alpha_star), nu1) + 2 * (1 - alpha) * (pt(pmax(q, 0)/(2*(1- alpha_star)), nu2) - 0.5)
  probs
}

#' @rdname astDist
#' @export
qast = function(p, mu=0, sigma=pi, alpha=0.5, nu1, nu2) {
  if(!is.numeric(p)) stop("p must be numeric")
  alpha_star = alpha*K(nu1) / (alpha * K(nu1)+(1-alpha)*K(nu2))
  quans = 2 * alpha_star * (mu + sigma * qt(pmin(p, alpha)/(2*alpha), nu1)) + 2 * (1 - alpha_star) * (mu + sigma * qt((pmax(p, alpha)+1-2*alpha)/(2*(1- alpha)), nu2))
  quans
}

#' @rdname astDist
#' @export
rast = function(n=1, mu=0, sigma=pi, alpha=0.5, nu1, nu2) {
  alpha_star = alpha*K(nu1) / (alpha * K(nu1)+(1-alpha)*K(nu2))
  u = runif(n)
  t1 = rt(n, nu1)
  t2 = rt(n, nu2)
  rvs = alpha_star * abs(t1) * (sign(u - alpha) - 1) +  (1 - alpha_star) * abs(t2) * (sign(u - alpha) + 1)
  rvs = mu + rvs * sigma * (alpha * K(nu1)+(1-alpha)*K(nu2))
  rvs
}

K = function(nu) {
  gamma(0.5*(nu+1)) / (sqrt(pi*nu)*gamma(0.5*nu))
}

check_bound = function(mu, sigma, alpha, nu1, nu2) {
  if(!is.numeric(mu)) stop("mu must be numeric")
  if(!is.numeric(sigma)) stop("mu must be numeric")
  if(!is.numeric(alpha)) stop("mu must be numeric")
  if(!is.numeric(nu1)) stop("mu must be numeric")
  if(!is.numeric(nu2)) stop("mu must be numeric")
  if(sigma <= 0) stop("sigma must be greater than 0")
  if(nu1 <= 0) stop("nu1 must be greater than 0")
  if(nu2 <= 0) stop("nu2 must be greater than 0")
  if(alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
}
