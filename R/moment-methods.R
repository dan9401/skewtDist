# Moment methods for AST & GAT distributions, while gat methods are not yet implemented
# we may also want to keep separate files for both distributions, doesn't seem necessary at the time
# and moment methods for ast & gat class without data, just for exploration uses
# may also want separate functions for mean, variance, sd, skewness & kurtosis
# authorized domain, here or in the plot file

moment_ast <- function(n, mu, sigma, alpha, nu1, nu2, method = c("numerical", "analytical")) {
  method <- match.arg(method)
  if (method == "analytical") {
    if (mu != 0) {
      stop("Analytical formula of moments cannot calculate with location parameters other than 0.")
    }
    # return value - analytical
    moment_ast_analytical(n, mu, sigma, alpha, nu1, nu2)
  } else {
    # return value - numerical
    moment_ast_numerical(n, mu, sigma, alpha, nu1, nu2)
  }
}

moment_ast_analytical <- function(n, sigma, alpha, nu1, nu2) {
  B <- alpha * K(nu1) + (1 - alpha) * K(nu2)
  alpha_star <- alpha * K(nu1)/B
  # return value
  alpha * (-2 * alpha_star * sigma * B)^n * moment_abs_t(nu1, n) +
    (1 - alpha) * (2 * (1 - alpha_star) * sigma * B)^n * moment_abs_t(nu2, n)
}

# needed change
moment_ast_numerical <- function(n, mu, sigma, alpha, nu1, nu2) {
  integrand <- function(x, n, mu, sigma, alpha, nu1, nu2) {
    x^n * dast(x, mu, sigma, alpha, nu1, nu2)
  }
  # return value
  # should use safeIntegrate instead, need testing
  integrate(integrand, -Inf, Inf, n, mu, sigma, alpha, nu1, nu2)
}
