#' @title Asymmetric Student t-distribution Functions
#'
#' @description Probablity density function (pdf), Cumulative distribution function (cdf), quantile function and random generation of AST
#' @param
#' @keywords
#' @examples
#' dast()
#' past()
#' qast()
#' rast()
#' @aliases past qast rast
#' @name ast

#' @rdname ast
#' @export
dast = function(x, mu, sigma, alpha, nu1, nu2) {
  K = function(nu) {
    gamma(0.5*(nu+1)) / (sqrt(pi*nu)*gamma(0.5*nu))
  }
  x1 = x[x <= mu]
  x2 = x[x > mu]
  dens = numeric(length(x))
  dens[x <= mu] = (1 + ((x1-mu)/(2*alpha*sigma*K(nu1)))^2/nu1)^(-0.5*(nu1+1)) / sigma
  dens[x > mu] = (1 + ((x2-mu)/(2*(1-alpha)*sigma*K(nu2)))^2/nu2)^(-0.5*(nu2+1)) / sigma
  dens
}

#' @rdname ast
#' @export
past = function(x, location, scale, skewness, dof1, dof2) {
  K = function(dof) {
    gamma(0.5*(dof+1)) / (sqrt(pi*dof)*gamma(0.5*dof))
  }
  alpha_star = skewness*K(dof1) / (skewness * K(dof1)+(1-skewness)*K(dof2))
  probs = 2 * skewness * pt(min(x, 0)/(2*alpha_star), dof1) + 2 * (1 - skewness) * (pt(max(x, 0)/(2*(1- alpha_star)), dof2) - 0.5)
  probs
}

#' @rdname ast
#' @export
qast = function(x, location, scale, skewness, dof1, dof2) {
  K = function(dof) {
    gamma(0.5*(dof+1)) / (sqrt(pi*dof)*gamma(0.5*dof))
  }
  alpha_star = skewness*K(dof1) / (skewness * K(dof1)+(1-skewness)*K(dof2))
  quans = 2 * alpha_star * qt(min(x, skewness)/(2*skewness), dof1) + 2 * (1 - alpha_star) * qt((max(x, skewness)+1-2*skewness)/(2*(1- skewness)), dof2)
  quans
}

#' @rdname ast
#' @export
rast = function(x, location, scale, skewness, dof1, dof2) {
  K = function(dof) {
    gamma(0.5*(dof+1)) / (sqrt(pi*dof)*gamma(0.5*dof))
  }
  alpha_star = skewness*K(dof1) / (skewness * K(dof1)+(1-skewness)*K(dof2))
  u = runif(x)
  t1 = rt(x, dof1)
  t2 = rt(x, dof2)
  rvs = alpha_star * abs(t1) * (sign(u - skewness) - 1) +  (1 - alpha_star) * abs(t2) * (sign(u - skewness) + 1)
  rvs
}
