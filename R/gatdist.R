#' @title Generalized Asymmetric t-distribution
#'
#' @description Probablity density function (pdf), Cumulative distribution function (cdf), quantile function and random generation of GAT distributions
#' @param x,q vector of quantiles
#' @param p vector of probablilities
#' @param n number of observations for random generation
#' @param mu location parameter
#' @param sigma scale parameter
#' @param r tail power asymmetry
#' @param c scale asymmetry
#' @param alpha how early tail behavior is apparent
#' @param nu degrees of freedom / tail parameter
#' @aliases dgat
#' @aliases pgat
#' @aliases rgat
#' @name gatDist

#' @rdname gatDist
#' @export
dgat = function(x, mu, sigma, alpha, r, c, nu) {
  if(!is.numeric(x)) stop("x must be numeric")
  g = (x-mu)/sigma + sqrt(1+((x-mu)/sigma)^2)
  A = alpha * (1 + r^2) / (r * sigma)
  B = (((c*g)^(alpha*r))+((c*g)^(-alpha/r)))^(-nu/alpha)/beta(nu/alpha/(1+r^2), r^2*nu/alpha/(1+r^2))
  C = (1+((x-mu)/sigma)^2)^(-0.5)
  d = A * B * C
  d
}

#' @rdname gatDist
#' @export
pgat = function(q, mu, sigma, alpha, r, c, nu){
  if(!is.numeric(x)) stop("x must be numeric")
  q = 1 / (1 + c^(-alpha*(1 + r^2)/r)*(((x-mu)/sigma)+sqrt(1+(x-mu)^2/sigma^2))^(-alpha*(1+r^2)/r))
  p = pbeta(q, nu/alpha/(1+r^2), r^2*nu/alpha/(1+r^2))
  p
}

#' @rdname gatDist
#' @export
qgat = function(p, mu, sigma, alpha, r, c, nu){
  if(!is.numeric(p)) stop("p must be numeric")
  if(p < 0 || p > 1) stop("p must be in (0,1)")
  q = 1 / (1 + c^(-alpha*(1 + r^2)/r)*(((x-mu)/sigma)+sqrt(1+(x-mu)^2/sigma^2))^(-alpha*(1+r^2)/r))
  p = pbeta(q, nu/alpha/(1+r^2), r^2*nu/alpha/(1+r^2))
  p
}

#' @rdname gatDist
#' @export
rgat = function(n, mu, sigma, alpha, r, c, nu) {
  if(n < 0) stop("x must be non-negative")
  a = nu/alpha/(1+r^2)
  b = nu*r^2/alpha/(1+r^2)
  q = rbeta(n, a, b)
  del = r / alpha / (1 + r^2)
  x = mu + 0.5*sigma*((q/(1-q))^del/c - c*(q/(1-q))^(-del))
  x
}
