#' @title Information Matrix Function of Asymmetric Student-t distribution
#'
#' @name gat-infomat
#' @aliases gatInfoMat
#'
#' @description Information matrix, asymptotic covariance and correlation matrix functions of Asymmetric Student-t distribution
#'
#' @param pars a vector of parameter values for an AST distribution
#' @param data a vector of numeric data used to calculate observed information matrix
#'
#' @details
#' The expected information matrix is calculated by the expectation of the outer product of score functions,
#' analytical formulas are provided in \emph{Zhu and Galbraith(2010)}.
#' The observed information matrix is calculated by the expectation of negative Hessian Matrix of the log-likelihood function.
#'
#' @examples
#' pars <- c(0.12, 0.6, 1.5, 1.2, 2, 5)
#' data <- rgat(1000, pars = pars)
#' 
#' round(gatInfoMat(pars, data), 4)
#' 
#' round(gatCov(pars, data), 4)
#' round(gatCor(pars, data), 4)
#' 
#' @importFrom stats cov2cor

#' @rdname gat-infomat
#' @export
gatInfoMat <- function(pars, data) {
  y <- data

  mu <- pars[1]
  phi <- pars[2]
  alpha <- pars[3]
  r <- pars[4]
  c <- pars[5]
  nu <- pars[6]
  T_ <- length(y)

  z <- (y - mu) / phi
  g <- z + sqrt(1 + z^2)
  p <- nu / (alpha * (1 + r^2))
  q <- p * r^2
  A <- (c*g)^(alpha*r) + (c*g)^(-alpha/r)
  dAdg <- alpha*r*(c*g)^(alpha*r)/g - alpha/r*(c*g)^(-alpha/r)/g
  dAdr <- (c*g)^(alpha*r)*alpha*log(c*g) + (c*g)^(-alpha/r)*alpha*log(c*g)/r^2
  dgdz <- 1 + z/sqrt(1+z^2)
  dzdmu <- -1/phi
  dzdphi <- -(y-mu)/phi^2
  dpdalpha <- -p/alpha
  dpdnu <- p/nu
  dpdr <- -2*r/(1+r^2)*p
  dqdalpha <- -q/alpha
  dqdnu <- q/nu
  dqdr <- 2/(r*(1+r^2))*q
  dbdp <- beta(p,q)*(digamma(p) - digamma(p+q))
  dbdq <- beta(p,q)*(digamma(q) - digamma(p+q))
  dbdnu <- dbdq * dqdnu + dbdp * dpdnu
  dbdalpha <- dbdq * dqdalpha + dbdp * dpdalpha
  dbdr <- dbdq * dqdr + dbdp * dpdr

  s_mu <- -nu/alpha/A*dAdg*dgdz*dzdmu - z/(1+z^2)*dzdmu
  s_phi<- -1/phi -nu/alpha/A*dAdg*dgdz*dzdphi - z/(1+z^2)*dzdphi
  s_alpha <-  1/alpha + nu/alpha^2*log(A) - nu/alpha/A*( r*(c*g)^(alpha*r)*log(c*g) - 1/r*(c*g)^(-alpha/r)*log(c*g) ) - 1/beta(p, q) * dbdalpha
  s_r <- 2*r/(1+r^2) - 1/r - nu/alpha/A*dAdr - 1/beta(p,q) * dbdr
  s_c <- -nu/alpha/A* (alpha*r*g^(alpha*r)*c^(alpha*r-1) - alpha/r*g^(-alpha/r)*c^(-alpha/r-1))
  s_nu <- -1/alpha*log(A) - 1/beta(p,q) * dbdnu

  e_mumu <- mean(s_mu * s_mu)
  e_muphi <- mean(s_mu * s_phi)
  e_mualpha <- mean(s_mu * s_alpha)
  e_mur <- mean(s_mu * s_r)
  e_muc <- mean(s_mu * s_c)
  e_munu <- mean(s_mu * s_nu)

  e_phiphi <- mean(s_phi * s_phi)
  e_phialpha <- mean(s_phi * s_alpha)
  e_phir <- mean(s_phi * s_r)
  e_phic <- mean(s_phi * s_c)
  e_phinu <- mean(s_phi * s_nu)

  e_alphaalpha <- mean(s_alpha * s_alpha)
  e_alphar <- mean(s_alpha * s_r)
  e_alphac <- mean(s_alpha * s_c)
  e_alphanu <- mean(s_alpha * s_nu)

  e_rr <- mean(s_r * s_r)
  e_rc <- mean(s_r * s_c)
  e_rnu <- mean(s_r * s_nu)
  e_cc <- mean(s_c * s_c)
  e_cnu <- mean(s_c * s_nu)
  e_nunu <- mean(s_nu * s_nu)

  infoMat <- matrix(c(e_mumu, e_muphi, e_mualpha, e_mur, e_muc, e_munu,
                      e_muphi, e_phiphi, e_phialpha, e_phir, e_phic, e_phinu,
                      e_mualpha, e_phialpha, e_alphaalpha, e_alphar, e_alphac, e_alphanu,
                      e_mur, e_phir, e_alphar, e_rr, e_rc, e_rnu,
                      e_muc, e_phic, e_alphac, e_rc, e_cc, e_cnu,
                      e_munu, e_phinu, e_alphanu, e_rnu, e_cnu, e_nunu),

                    nrow = 6, ncol = 6)

  infoMat
}

#' @rdname gat-infomat
#' @export
gatCov <- function(pars, data) {
  solve(gatInfoMat(pars, data))
}

#' @rdname gat-infomat
#' @export
gatCor <- function(pars, data) {
  cov2cor(gatCov(pars, data))
}
