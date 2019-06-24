#' @title Asymmetric Student t-distribution
#'
#' @description Probablity density function (pdf), Cumulative distribution function (cdf), quantile function and random generation of AST distributions
#' @param x,q vector of quantiles
#' @param p vector of probablilities
#' @param n number of observations for random generation
#' @param mu location parameter
#' @param sigma scale parameter
#' @param alpha skewness parameter
#' @param nu1 degrees of freedom / tail parameter 1
#' @param nu2 degrees of freedom / tail parameter 2
#'
#' @name infoMat_ast
#' @examples

#' @rdname infoMat
#' @export
infoMat_ast <- function(fit, method = c("expected", "observed")) {
  y <- fit$data
  pars <- fit$fitted_pars
  method <- match.arg(method, c("expected", "observed"))
  mu <- pars["mu"]
  sigma <- pars["sigma"]
  alpha <- pars["alpha"]
  nu1 <- pars["nu1"]
  nu2 <- pars["nu2"]
  T_ <- length(y)
  y1 <- y[y <= mu]
  y2 <- y[y > mu]

  if (method == "expected") {
    S_mu2 <- 1/(4*sigma^2) * ( (nu1+1)/(alpha*(nu1+3))/K(nu1)^2 + (nu2+1)/((1-alpha)*(nu2+3))/K(nu2)^2 )
    S_sigma2 <- 2/sigma^2 * ( alpha*nu1/(nu1+3) + (1-alpha)*nu2/(nu2+3) )
    S_alpha2 <- 3 * ( (nu1+1)/(alpha*(nu1+3)) + (nu2+1)/((1-alpha)*(nu2+3))  )
    S_nu12 <- alpha/2 * ( nu1/(nu1+3)*D(nu1)^2 - 2/(nu1+1)*D(nu1) - Dprime(nu1) )
    S_nu22 <- (1-alpha)/2 * ( nu2/(nu2+3)*D(nu2)^2 - 2/(nu2+1)*D(nu2) - Dprime(nu2) )

    S_musigma <- 2/sigma^2 * ( -(nu1+1)/(nu1+3) + (nu2+1)/(nu2+3) )
    S_mualpha <- -2/sigma * ( (nu1+1)/(alpha*(nu1+3)) + (nu2+1)/((1-alpha)*(nu2+3)) )
    S_munu1 <- 1/sigma * ( 1/(nu1+1) - (nu1+1)*D(nu1)/(nu1+3) )
    S_munu2 <- - 1/sigma * ( 1/(nu2+1) - (nu2+1)*D(nu2)/(nu2+3) )

    S_sigmaalpha <- 2/sigma * ( nu1/(nu1+3) - nu2/(nu2+3) )
    S_sigmanu1 <- alpha/sigma * ( -1/(nu1+1) + nu1*D(nu1)/(nu1+3) )
    S_sigmanu2 <- (1-alpha)/sigma * ( -1/(nu2+1) + nu2*D(nu2)/(nu2+3) )

    S_alphanu1 <- -1/(nu1+1) + nu1*D(nu1)/(nu1+3)
    S_alphanu2 <- 1/(nu2+1) - nu2*D(nu2)/(nu2+3)
    S_nu1nu2 <- 0
  } else {
    ###
    S_mu2 <- 1/nu1 * ((nu1+1)/(2*alpha*sigma*K(nu1)))^2 * mean(1/L(pars, y1) - 1/L(pars, y1)^2) +
      1/nu2 * ((nu2+1)/(2*(1-alpha)*sigma*K(nu2)))^2 * mean(1/R(pars, y2) - 1/R(pars, y2)^2)
    S_sigma2 <- -1/sigma^2 + ((nu1+1)/sigma)^2 * mean((1 - 1/L(pars, y1))^2) +
      ((nu2+1)/sigma)^2 * mean((1 - 1/R(pars, y2))^2)
    S_alpha2 <- ((nu1+1)/alpha)^2 * mean((1 - 1/L(pars, y1))^2) +
      ((nu2+1)/(1-alpha))^2 * mean((1 - 1/R(pars, y2))^2)
    S_nu12 <- mean( log(L(pars, y1))^2/4 + ((nu1+1)/2)^2*D(nu1)^2*(1 - 1/L(pars, y1))^2 ) -
      (nu1+1)/2*D(nu1) * mean( (1 - 1/L(pars, y1))*log(L(pars, y1)) )
    S_nu22 <- mean( log(R(pars, y2))^2/4 + ((nu2+1)/2)^2*D(nu2)^2*(1 - 1/R(pars, y2))^2 ) -
      (nu2+1)/2*D(nu2) * mean( (1 - 1/R(pars, y2))*log(R(pars, y2)) )

    S_musigma <- (nu1+1)^2/(2*sigma) * mean( (1/L(pars, y1) - 1/L(pars, y1)^2) / nu1*2*(y1-mu)/(2*alpha*K(nu1))^2) +
      (nu2+1)^2/(2*sigma) * mean( (1/R(pars, y2) - 1/R(pars, y2)^2) / nu2*2*(y2-mu)/(2*(1-alpha)*K(nu2))^2)
    S_mualpha <- (nu1+1)^2/(2*alpha) * mean( (1/L(pars, y1) - 1/L(pars, y1)^2) / nu1*2*(y1-mu)/(2*alpha*sigma*K(nu1))^2) -
      (nu2+1)^2/(2*(1-alpha)) * mean( (1/R(pars, y2) - 1/R(pars, y2)^2) / nu2*2*(y2-mu)/(2*(1-alpha)*sigma*K(nu2))^2)
    # mu nu1 may not be accurate
    S_munu1 <- -(nu1+1)/4 * mean( 1/L(pars, y1)/nu1 * 2*(y1-mu)/(2*alpha*sigma*K(nu1))^2 * log(L(pars, y1)) ) +
      ((nu1+1)/2)^2*D(nu1) * mean( 1/L(pars, y1)/nu1 * 2*(y1-mu)/(2*alpha*sigma*K(nu1))^2 * (1 - 1/L(pars, y1)) )
    S_munu2 <- (nu2+1)/4 * mean( 1/R(pars, y2)/nu2 * 2*(y2-mu)/(2*(1-alpha)*sigma*K(nu2))^2 * log(R(pars, y2)) ) +
      ((nu2+1)/2)^2*D(nu2) * mean( 1/R(pars, y2)/nu2 * 2*(y2-mu)/(2*(1-alpha)*sigma*K(nu2))^2 * (1 - 1/R(pars, y2)))

    S_sigmaalpha <- (nu1+1)^2/(sigma*alpha) * mean((1 - 1/L(pars, y1))^2) -
      (nu2+1)^2/(sigma*(1-alpha)) * mean((1 - 1/R(pars, y2))^2)
    S_sigmanu1 <- -(nu1+1)/(2*alpha) * mean( (1 - 1/L(pars, y1)) * log(L(pars, y1)) ) +
      (nu1+1)^2/(2*alpha) * D(nu1) * mean( (1 - 1/L(pars, y1))^2 )
    S_sigmanu2 <- -(nu2+1)/(2*(1-alpha)) * mean( (1 - 1/R(pars, y2)) * log(L(pars, y2)) ) +
      (nu2+1)^2/(2*(1-alpha)) * D(nu2) * mean( (1 - 1/R(pars, y2))^2 )

    S_alphanu1 <- -(nu1+1)/(2*alpha) * mean( (1 - 1/L(pars, y1))*log(L(pars, y1)) ) +
      (nu1+1)^2/(2*alpha) * D(nu1) * mean((1 - 1/L(pars, y1))^2)
    S_alphanu2 <- (nu2+1)/(2*(1-alpha)) * mean( (1 - 1/R(pars, y2))*log(R(pars, y2)) ) -
      (nu2+1)^2/(2*(1-alpha)) * D(nu2) * mean((1 - 1/R(pars, y2))^2)
    S_nu1nu2 <- 0


  }

  infoMat <- matrix(c(S_mu2, S_musigma, S_mualpha, S_munu1, S_munu2,
                      S_musigma, S_sigma2, S_sigmaalpha, S_sigmanu1, S_sigmanu2,
                      S_mualpha, S_sigmaalpha, S_alpha2, S_alphanu1, S_alphanu2,
                      S_munu1, S_sigmanu1, S_alphanu1, S_nu12, S_nu1nu2,
                      S_munu2, S_sigmanu2, S_alphanu2, S_nu1nu2, S_nu22),
                    nrow = 5, ncol = 5)
  rownames(infoMat) = colnames(infoMat) = c("mu", "sigma", "alpha", "nu1", "nu2")
  infoMat
}

infoMat_t <- function(sigma, nu) {
  I11 <- 1/4*(trigamma((nu+1)/2) - trigamma(nu/2)) - 1/nu*(1/(nu+1) - 1/(2*(nu+3)) )
  I12 <- 1/sigma*(1/(nu+3)-1/(nu+1))
  I22 <- 2/sigma^2*nu/(nu+3)
  I13 <- I23 <- 0
  I33 <- 1/sigma^2*(1-2/(nu+3))
  infoMat <- matrix(c(I11, I12, I13,
                      I12, I22, I23,
                      I13, I23, I33),
                    nrow = 3, ncol = 3)
  rownames(infoMat) = colnames(infoMat) = c("nu", "sigma", "mu")
  infoMat
}
