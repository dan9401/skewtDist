gatfit <- function(data, start_pars = c(), fixed_pars = c(), solver = c("nloptr", "Rsolnp"), solver_control) {




}

llgat <- function(pars, arglist) {
  y <- arglist$data
  fixed_pars <- arglist$fixed_pars

  est_idx <- which(is.na(fixed_pars))
  fix_idx <- which(!(is.na(fixed_pars)))
  all_pars <- c()
  all_pars[est_idx] <- pars
  all_pars[fix_idx] <- fixed_pars[fix_idx]

  mu <- all_pars[1]
  phi <- all_pars[2]
  alpha <- all_pars[3]
  r <- all_pars[4]
  c <- all_pars[5]
  nu <- all_pars[6]
  T_ <- length(y)

  z <- (y - mu) / phi
  g <- z + sqrt(1 + z^2)
  p <- nu / (alpha * (1 + r^2))
  q <- p * r^2

  logl <- T_ * log( alpha * (1 + r^2) / (r * phi) ) - sum( nu / alpha * log( (cg)^(alpha*r) + cg^(-alpha/r) ) ) - T_ * log(beta(p, q)) - sum( 1/2 * log(1 + z^2) )
  -logl
}

llgat_grad <- function(pars, arglist) {
  y <- arglist$data
  fixed_pars <- arglist$fixed_pars

  est_idx <- which(is.na(fixed_pars))
  fix_idx <- which(!(is.na(fixed_pars)))
  all_pars <- c()
  all_pars[est_idx] <- pars
  all_pars[fix_idx] <- fixed_pars[fix_idx]

  mu <- all_pars[1]
  phi <- all_pars[2]
  alpha <- all_pars[3]
  r <- all_pars[4]
  c <- all_pars[5]
  nu <- all_pars[6]
  T_ <- length(y)

  z <- (y - mu) / phi
  g <- z + sqrt(1 + z^2)
  p <- nu / (alpha * (1 + r^2))
  q <- p * r^2

  g_mu <- sum((nu1 + 1) / L(all_pars, y1) / nu1 * (y1 - mu) / (2 * alpha * phi * K(nu1))^2) +
    sum((nu2 + 1) / R(all_pars, y2) / nu2 * (y2 - mu) / (2 * (1 - alpha) * phi * K(nu2))^2)
  g_phi<- -T_ / phi + (nu1 + 1) / phi * sum(1 - 1 / L(all_pars, y1)) + (nu2 + 1) / phi * sum(1 - 1 / R(all_pars, y2))
  g_alpha <- (nu1 + 1) / alpha * sum(1 - 1 / L(all_pars, y1)) - (nu2 + 1) / (1 - alpha) * sum(1 - 1 / R(all_pars, y2))
  g_r <- sum(- log(L(all_pars, y1)) / 2 + (nu1 + 1) / 2 * D(nu1) * (L(all_pars, y1) - 1) / L(all_pars, y1))
  g_c <- sum(- log(R(all_pars, y2)) / 2 + (nu2 + 1) / 2 * D(nu2) * (R(all_pars, y2) - 1) / R(all_pars, y2))
  g_nu <- sum(- log(R(all_pars, y2)) / 2 + (nu2 + 1) / 2 * D(nu2) * (R(all_pars, y2) - 1) / R(all_pars, y2))
  gradient <- -c(mu = g_mu, phi = g_phi, alpha = g_alpha, g_r = g_r, g_c = g_c, nu = g_nu)

  return(gradient[est_idx])
}
