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
#' @param x,q vector of quantiles
#' @param p vector of probablilities
#' @param n number of observations for random generation
#' @param mu location parameter
#' @param s scale parameter, \eqn{s > 0}
#' @param alpha skewness parameter, \eqn{0 < alpha < 1}
#' @param nu degrees of freedom / tail parameter for the both tails, \eqn{ nu > 0}
#' @param pars a vector that contains mu, s, alpha, nu1, nu2, if pars is specified, mu, s, alpha, nu1, nu2 should not be specified
#' @param method method used to calculate the moment(s), one of 'analytical' and 'numerical'
#' @param type type of kurtosis calculated, one of 'excess' and 'regular'
#' 
#' @param data a univariate data object to be fitted
#' @param start_pars a named numeric vector of starting parameters for the optimization algorithm, not all parameters are needed
#' @param fixed_pars a named numeric vector of parameters to be kept fixed during the optimization routine, not all parameters are needed
#' @param solver solver used for MLE, one of 'nlminb', 'nloptr', 'Rsolnp', default is 'nlminb'
#' @param solver_control list of control arguments passed to the solver
#'
#' @details
#' The SST functions are wrappers of their AST equivalent.
#'
#' For most of the functions, they are implemented by simply setting nu1, nu2(AST) = nu(SST).
#'
#' For \code{sstFit}, it is implemented by setting symmetric = TRUE.

#' @rdname sst-functions
#' @export
dsst <- function(x, mu = 0, s = 1, alpha = 0.5, nu = Inf, pars = NULL) {
  dast(x, mu, s, alpha, nu, nu, pars)
}

#' @rdname sst-functions
#' @export
psst <- function(q, mu = 0, s = 1, alpha = 0.5, nu = Inf, pars = NULL) {
  past(q, mu, s, alpha, nu, nu, pars)
}

#' @rdname sst-functions
#' @export
qsst <- function(p, mu = 0, s = 1, alpha = 0.5, nu = Inf, pars = NULL) {
  qast(p, mu, s, alpha, nu, nu, pars)
}

#' @rdname sst-functions
#' @export
rsst <- function(n, mu = 0, s = 1, alpha = 0.5, nu = Inf, pars = NULL) {
  rast(n, mu, s, alpha, nu, nu, pars)
}

#' @rdname sst-functions
#' @export
sstMoments <- function(mu = 0, s = 1, alpha = 0.5, nu = Inf, pars = NULL, method = c("analytical", "numerical"), type = c("excess", "regular")) {
  astMoments(mu, s, alpha, nu, nu, pars, method, type)
}

#' @rdname sst-functions
#' @export
sstRawMoment <- function(n, mu = 0, s = 1, alpha = 0.5, nu = Inf, pars = NULL, method = c("analytical", "numerical")) {
  astRawMoment(n, mu, s, alpha, nu, nu, pars, method)
}

#' @rdname sst-functions
#' @export
sstCentralMoment <- function(n, mu = 0, s = 1, alpha = 0.5, nu = Inf, pars = NULL, method = c("analytical", "numerical")) {
  astCentralMoment(n, mu, s, alpha, nu, nu, pars, method)
}

#' @rdname sst-functions
#' @export
sstInfoMat <- function(pars) {
  mu <- pars[1]
  s <- pars[2]
  alpha <- pars[3]
  nu <- pars[4]

  S_mu2 <- (nu + 1) / ((4*s^2) * alpha * (1 - alpha) * (nu + 3) * K(nu)^2 )
  S_s2 <- 2 * nu / (s^2 * (nu+3) )
  S_alpha2 <- 3 * (nu+1) / (alpha * (1-alpha) * (nu+3) )
  S_nu2 <- 1/2 * ( nu/(nu+3)*D(nu)^2 - 2/(nu+1)*D(nu) - Dprime(nu) )

  S_mualpha <- -2/s * (nu+1) / (alpha * (1-alpha) * (nu+3) )
  S_s_nu <- 1/s * ( -1/(nu+1) + nu*D(nu)/(nu+3) )

  S_mu_s <- 0
  S_munu <- 0
  S_s_alpha <- 0
  S_alphanu <- 0

  infoMat <- matrix(c(S_mu2, S_mu_s, S_mualpha, S_munu,
                      S_mu_s, S_s2, S_s_alpha, S_s_nu,
                      S_mualpha, S_s_alpha, S_alpha2, S_alphanu,
                      S_munu, S_s_nu, S_alphanu, S_nu2),
                    nrow = 4, ncol = 4)
  rownames(infoMat) = colnames(infoMat) = c("mu", "s", "alpha", "nu")
  infoMat
}

#' @rdname sst-functions
#' @export
sstMLE <- function(data, start_pars = c(), fixed_pars = c(), solver = c("nlminb", "nloptr", "Rsolnp"), solver_control = list()) {
  astMLE(data, start_pars, fixed_pars, solver, solver_control, symmetric = TRUE)
}


