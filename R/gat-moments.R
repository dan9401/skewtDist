mean_gat <- function(mu, phi, alpha, r, c, nu) {
  A <- nu/( alpha*(1+r^2) )
  B <- A * r^2
  delta <- A / nu * r
  mu + phi * (beta(A + delta, B - delta)/c - c*beta(A - delta, B + delta)) /
                (2*beta(A, B))
}

moment_central_gat <- function(n, mu, phi, alpha, r, c, nu) {
  A <- nu/( alpha*(1+r^2) )
  B <- A * r^2
  delta <- A / nu * r
  m <- 0:n
  (phi/2)^n / beta(A, B) * sum( (-1)^m * choose(n, m) * c^(n - 2*m) *
                                  beta(A - (n-2*m)*delta, B + (n-2*m)*delta))
}

var_gat <- function(mu, phi, alpha, r, c, nu) {
  A <- nu/( alpha*(1+r^2) )
  B <- A * r^2
  delta <- A / nu * r
  phi^2 / (4*beta(A, B)) * (c^(-2)*beta(A + 2*delta, B - 2*delta) + c^2*beta(A - 2*delta, B + 2* delta)) -
    phi^2/2
}

# var_gat = moment_central_gat but not var(data)
vg <- function(mu, phi, alpha, r, c, nu) {
  m <- mean_gat(mu, phi, alpha, r, c, nu)
  integrand <- function(x) {
    (x - m)^2 * dgat(x, mu, phi, alpha, r, c, nu)
  }
  safeIntegrate(integrand, -Inf, Inf)$value
}

sg <- function(mu, phi, alpha, r, c, nu) {
  m <- mean_gat(mu, phi, alpha, r, c, nu)
  sd <- sqrt(vg(mu, phi, alpha, r, c, nu))
  integrand <- function(x) {
    ((x - m)/sd)^3 * dgat(x, mu, phi, alpha, r, c, nu)
  }
  safeIntegrate(integrand, -Inf, Inf)$value
}
