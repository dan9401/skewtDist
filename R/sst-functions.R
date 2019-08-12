#' @title Functions for Symmetric Student-t distribution
#'
#' @name sst-functions
#' @aliases sst
#' @aliases dsst
#' @aliases psst
#' @aliases qsst
#' @aliases rsst
#' @aliases sstMoment
#' @aliases sstMoments
#' @aliases sstRawMoments
#' @aliases sstCentralMoments
#' @aliases sstMean
#' @aliases sstVar
#' @aliases sstSD
#' @aliases sstSkew
#' @aliases sstKurt
#' @aliases sstInfoMat
#' @aliases sstFit
#'
#' @details
#' The SST functions are wrappers of their AST equivalent.
#'
#' For most of the functions, they are implemented by simply setting nu1, nu2(AST) = nu(SST).
#'
#' For \code{sstFit}, it is implemented by setting symmetric = TRUE.

#' @rdname sst-functions
#' @export
dsst <- function(x, mu = 0, sigma = 1, alpha = 0.5, nu = Inf, pars = NULL) {
  dast(x, mu, sigma, alpha, nu, nu, pars)
}

#' @rdname sst-functions
#' @export
psst <- function(q, mu = 0, sigma = 1, alpha = 0.5, nu = Inf, pars = NULL) {
  past(q, mu, sigma, alpha, nu, nu, pars)
}

#' @rdname sst-functions
#' @export
qsst <- function(p, mu = 0, sigma = 1, alpha = 0.5, nu = Inf, pars = NULL) {
  qast(p, mu, sigma, alpha, nu, nu, pars)
}

#' @rdname sst-functions
#' @export
rsst <- function(n, mu = 0, sigma = 1, alpha = 0.5, nu = Inf, pars = NULL) {
  rast(n, mu, sigma, alpha, nu, nu, pars)
}

#' @rdname sst-functions
#' @export
sstMoment <- function(moment, mu = 0, sigma = 1, alpha = 0.5, nu = Inf, pars = NULL, method = c("analytical", "numerical"), type = c("excess", "regular")) {
  astMoment(moment, mu, sigma, alpha, nu, nu, pars, method, type)
}

#' @rdname sst-functions
#' @export
sstMoments <- function(mu = 0, sigma = 1, alpha = 0.5, nu = Inf, pars = NULL, method = c("analytical", "numerical"), type = c("excess", "regular")) {
  astMoments(mu, sigma, alpha, nu, nu, pars, method, type)
}

#' @rdname sst-functions
#' @export
sstRawMoment <- function(n, mu = 0, sigma = 1, alpha = 0.5, nu = Inf, pars = NULL, method = c("analytical", "numerical")) {
  astRawMoment(n, mu, sigma, alpha, nu, nu, pars, method)
}

#' @rdname sst-functions
#' @export
sstCentralMoment <- function(n, mu = 0, sigma = 1, alpha = 0.5, nu = Inf, pars = NULL, method = c("analytical", "numerical")) {
  astCentralMoment(n, mu, sigma, alpha, nu, nu, pars, method)
}

#' @rdname sst-functions
#' @export
sstMean <- function(mu = 0, sigma = 1, alpha = 0.5, nu = Inf, pars = NULL, method = c("analytical", "numerical")) {
  astMean(mu, sigma, alpha, nu, nu, pars, method)
}

#' @rdname sst-functions
#' @export
sstSD <- function(mu = 0, sigma = 1, alpha = 0.5, nu = Inf, pars = NULL, method = c("analytical", "numerical")) {
  astSD(mu, sigma, alpha, nu, nu, pars, method)
}

#' @rdname sst-functions
#' @export
sstVar <- function(mu = 0, sigma = 1, alpha = 0.5, nu = Inf, pars = NULL, method = c("analytical", "numerical")) {
  astVar(mu, sigma, alpha, nu, nu, pars, method)
}

#' @rdname sst-functions
#' @export
sstSkew <- function(mu = 0, sigma = 1, alpha = 0.5, nu = Inf, pars = NULL, method = c("analytical", "numerical")) {
  astKurt(mu, sigma, alpha, nu, nu, pars, method)
}

#' @rdname sst-functions
#' @export
sstKurt <- function(mu = 0, sigma = 1, alpha = 0.5, nu = Inf, pars = NULL, method = c("analytical", "numerical"), type = c("excess", "regular")) {
  astKurt(mu, sigma, alpha, nu, nu, pars, method)
}

#' @rdname sst-functions
#' @export
sstInfoMat <- function(pars) {
  mu <- pars[1]
  sigma <- pars[2]
  alpha <- pars[3]
  nu <- pars[4]

  S_mu2 <- (nu + 1) / ((4*sigma^2) * alpha * (1 - alpha) * (nu + 3) * K(nu)^2 )
  S_sigma2 <- 2 * nu / (sigma^2 * (nu+3) )
  S_alpha2 <- 3 * (nu+1) / (alpha * (1-alpha) * (nu+3) )
  S_nu2 <- 1/2 * ( nu/(nu+3)*D(nu)^2 - 2/(nu+1)*D(nu) - Dprime(nu) )

  S_mualpha <- -2/sigma * (nu+1) / (alpha * (1-alpha) * (nu+3) )
  S_sigmanu <- 1/sigma * ( -1/(nu+1) + nu*D(nu)/(nu+3) )

  S_musigma <- 0
  S_munu <- 0
  S_sigmaalpha <- 0
  S_alphanu <- 0

  infoMat <- matrix(c(S_mu2, S_musigma, S_mualpha, S_munu,
                      S_musigma, S_sigma2, S_sigmaalpha, S_sigmanu,
                      S_mualpha, S_sigmaalpha, S_alpha2, S_alphanu,
                      S_munu, S_sigmanu, S_alphanu, S_nu2),
                    nrow = 4, ncol = 4)
  rownames(infoMat) = colnames(infoMat) = c("mu", "sigma", "alpha", "nu")
  infoMat
}

#' @rdname sst-functions
#' @export
sstFit <- function(data, start_pars = c(), fixed_pars = c(), solver = c("nlminb", "nloptr", "Rsolnp"), solver_control = list()) {
  astFit(data, start_pars, fixed_pars, solver, solver_control, symmetric = TRUE)
}


