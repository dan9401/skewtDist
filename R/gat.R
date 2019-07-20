
gat <- function(mu, phi, alpha, r, c, nu) {
  structure(c(mu = mu, phi = phi, alpha = alpha, r = r, c = c, nu = nu), class = "gat")
}

dgat.gat <- function(x, gat) {
  mu <- gat[1]
  phi <- gat[2]
  alpha <- gat[3]
  r <- gat[4]
  c <- gat[5]
  nu <- gat[6]
  d <- dgat(x, mu, phi, alpha, r, c, nu)
  names(d) <- NULL
  d
}

pgat.gat <- function(gat, x) {
  mu <- gat[1]
  phi <- gat[2]
  alpha <- gat[3]
  r <- gat[4]
  c <- gat[5]
  nu <- gat[6]
  p <- pgat(x, mu, phi, alpha, r, c, nu)
  names(p) <- NULL
  p
}

qgat.gat <- function(gat, x) {
  mu <- gat[1]
  phi <- gat[2]
  alpha <- gat[3]
  r <- gat[4]
  c <- gat[5]
  nu <- gat[6]
  q <- qgat(x, mu, phi, alpha, r, c, nu)
  names(q) <- NULL
  q
}

rgat.gat <- function(gat, x) {
  mu <- gat[1]
  phi <- gat[2]
  alpha <- gat[3]
  r <- gat[4]
  c <- gat[5]
  nu <- gat[6]
  r <- rgat(x, mu, phi, alpha, r, c, nu)
  names(r) <- NULL
  r
}

dgat <- function(x, class, ...) {
  UseMethod("dgat", class)
}

pgat <- function(x, class...) {
  UseMethod("pgat", class)
}

qgat <- function(x, class...) {
  UseMethod("qgat", class)
}

rgat <- function(x, class...) {
  UseMethod("rgat", class)
}


