#' @title Generalized Asymmetric t-distribution
#'
#' @description Probablity density function (pdf), Cumulative distribution function (cdf), quantile function and random generation of GAT distributions
#' @param xvector of quantiles
#' @param n number of observations for random generation
#' @param mu location parameter
#' @param sigma scale parameter
#' @param r tail power asymmetry
#' @param c scale asymmetry
#' @param alpha how early tail behavior is apparent
#' @aliases dgat
#' @aliases pgat
#' @aliases rgat
#' @name gatDist

#' @rdname gatDist
#' @export
dgat = function(x, mu, sigma, dof, skewness, r, c) {
  g = (x-mu)/sigma + sqrt(1+((x-mu)/sigma)^2)
  A = skewness * (1 + r^2) / (r * sigma)
  B = (((c*g)^(skewness*r))+((c*g)^(-skewness/r)))^(-dof/skewness)/beta(dof/skewness/(1+r^2), r^2*dof/skewness/(1+r^2))
  C = (1+((x-mu)/sigma)^2)^(-0.5)
  A * B * C
}

#' @rdname gatDist
#' @export
pgat = function(x, mu, sigma, dof, skewness, r, c){
  q = 1 / (1 + c^(-skewness*(1 + r^2)/r)*(((x-mu)/sigma)+sqrt(1+(x-mu)^2/sigma^2))^(-skewness*(1+r^2)/r))
  pbeta(q, dof/skewness/(1+r^2), r^2*dof/skewness/(1+r^2))
}

#' @rdname gatDist
#' @export
rgat = function(n, mu, sigma, dof, skewness, r, c) {
  a = dof/skewness/(1+r^2)
  b = dof*r^2/skewness/(1+r^2)
  q = rbeta(n, a, b)
  del = r / skewness / (1 + r^2)
  x = mu + 0.5*sigma*((q/(1-q))^del/c - c*(q/(1-q))^(-del))
  x
}
