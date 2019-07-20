
ast <- function(mu, sigma, alpha, nu1, nu2) {
  structure(c(mu = mu, sigma = sigma, alpha = alpha, nu1 = nu1, nu2 = nu2), class = "ast")
}

dast.ast <- function(x, ast) {
  mu <- ast[1]
  sigma <- ast[2]
  alpha <- ast[3]
  nu <- ast[4]
  nu <- ast[5]
  d <- dast(x, mu, sigma, alpha, nu1, nu2)
  names(d) <- NULL
  d
}

past.ast <- function(x, ast) {
  mu <- ast[1]
  sigma <- ast[2]
  alpha <- ast[3]
  nu <- ast[4]
  nu <- ast[5]
  p <- past(x, mu, sigma, alpha, nu1, nu2)
  names(p) <- NULL
  p
}

qast.ast <- function(x, ast) {
  mu <- ast[1]
  sigma <- ast[2]
  alpha <- ast[3]
  nu <- ast[4]
  nu <- ast[5]
  q <- qast(x, mu, sigma, alpha, nu1, nu2)
  names(q) <- NULL
  q
}

rast.ast <- function(x, ast) {
  mu <- ast[1]
  sigma <- ast[2]
  alpha <- ast[3]
  nu <- ast[4]
  nu <- ast[5]
  r <- rast(x, mu, sigma, alpha, nu1, nu2)
  names(r) <- NULL
  r
}

dast <- function(x, class, ...) {
  UseMethod("dast", class)
}

past <- function(x, class...) {
  UseMethod("past", class)
}

qast <- function(x, class...) {
  UseMethod("qast", class)
}

rast <- function(x, class...) {
  UseMethod("rast", class)
}

moment.ast <- function() {}
