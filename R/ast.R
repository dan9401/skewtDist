
#' @export
ast <- function(mu, sigma, alpha, nu1, nu2) {
  structure(c(mu = mu, sigma = sigma, alpha = alpha, nu1 = nu1, nu2 = nu2), class = "ast")
}

#' @export
dast.ast <- function(x, ast) {
  mu <- ast[1]
  sigma <- ast[2]
  alpha <- ast[3]
  nu1 <- ast[4]
  nu2 <- ast[5]
  d <- dast(x, mu, sigma, alpha, nu1, nu2)
  names(d) <- NULL
  d
}

#' @export
past.ast <- function(x, ast) {
  mu <- ast[1]
  sigma <- ast[2]
  alpha <- ast[3]
  nu1 <- ast[4]
  nu2 <- ast[5]
  p <- past(x, mu, sigma, alpha, nu1, nu2)
  names(p) <- NULL
  p
}

#' @export
qast.ast <- function(x, ast) {
  mu <- ast[1]
  sigma <- ast[2]
  alpha <- ast[3]
  nu1 <- ast[4]
  nu2 <- ast[5]
  q <- qast(x, mu, sigma, alpha, nu1, nu2)
  names(q) <- NULL
  q
}

#' @export
rast.ast <- function(x, ast) {
  mu <- ast[1]
  sigma <- ast[2]
  alpha <- ast[3]
  nu1 <- ast[4]
  nu2 <- ast[5]
  r <- rast(x, mu, sigma, alpha, nu1, nu2)
  names(r) <- NULL
  r
}

#' @export
dast <- function(x, class, ...) {
  UseMethod("dast", class)
}

#' @export
past <- function(x, class, ...) {
  UseMethod("past", class)
}

#' @export
qast <- function(x, class, ...) {
  UseMethod("qast", class)
}

#' @export
rast <- function(x, class, ...) {
  UseMethod("rast", class)
}

moment.ast <- function() {}




