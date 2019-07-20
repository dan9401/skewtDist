# Moment methods for AST & GAT distributions, while gat methods are not yet implemented
# we may also want to keep separate files for both distributions, doesn't seem necessary at the time
# and moment methods for ast & gat class without data, just for exploration uses
# may also want separate functions for mean, variance, sd, skewness & kurtosis
# authorized domain, here or in the plot file



mean_ast <- function(mu, sigma, alpha, nu1, nu2) {
  integrand <- function(x) {
    x * dast(x, mu, sigma, alpha, nu1, nu2)
  }
  safeIntegrate(integrand, -Inf, Inf)$value
}

variance_ast <- function(mu, sigma, alpha, nu1, nu2) {
  mean <- mean_ast(mu, sigma, alpha, nu1, nu2)
  integrand <- function(x) {
    (x - mean)^2 * dast(x, mu, sigma, alpha, nu1, nu2)
  }
  safeIntegrate(integrand, -Inf, Inf)$value
}

sd_ast <- function(mu, sigma, alpha, nu1, nu2) {
  sqrt(variance_ast(mu, sigma, alpha, nu1, nu2))
}

skewness_ast <- function(mu, sigma, alpha, nu1, nu2) {
  mean <- mean_ast(mu, sigma, alpha, nu1, nu2)
  sd <- sd_ast(mu, sigma, alpha, nu1, nu2)
  integrand <- function(x) {
    ((x - mean) / sd)^3 * dast(x, mu, sigma, alpha, nu1, nu2)
  }
  safeIntegrate(integrand, -Inf, Inf)$value
}

kurtosis_ast <- function(mu, sigma, alpha, nu1, nu2) {
  mean <- mean_ast(mu, sigma, alpha, nu1, nu2)
  sd <- sd_ast(mu, sigma, alpha, nu1, nu2)
  integrand <- function(x) {
    ((x - mean) / sd)^4 * dast(x, mu, sigma, alpha, nu1, nu2)
  }
  safeIntegrate(integrand, -Inf, Inf)$value
}

moment_ast_analytical <- function(n, sigma, alpha, nu1, nu2) {
  B <- alpha * K(nu1) + (1 - alpha) * K(nu2)
  alpha_star <- alpha * K(nu1)/B
  # return value
  alpha * (-2 * alpha_star * sigma * B)^n * moment_abs_t(nu1, n) +
    (1 - alpha) * (2 * (1 - alpha_star) * sigma * B)^n * moment_abs_t(nu2, n)
}

moment_abs_t <- function(nu, n) {
  # -1 < n < nu
  sqrt(nu^n / pi) * gamma( (n+1)/2 ) * gamma( (nu-n)/2 ) / gamma(nu/2)
}

# needed change
moment_ast_numerical <- function(n, mu, sigma, alpha, nu1, nu2) {
  integrand <- function(x) {
    x^n * dast(x, mu, sigma, alpha, nu1, nu2)
  }
  safeIntegrate(integrand, -Inf, Inf)$value
}

