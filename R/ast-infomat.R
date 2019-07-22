#' @title Information Matrix for AST and GAT distributions
#'
#' @description Information Matrix function for AST and GAT distributions, function for GAT distributions have not yet been developed
#' we may want information matrix for symmetric t distributions, or functions for ast & gat objects without data, for users to explore the qualities of distribuitons
#' we may also want to keep separate files for both distributions, doesn't seem necessary at the time
#' and infoMat functions on fit and dist objects
#'
#' @param pars vector of parameter values for an AST distribution (or an GAT distribution)
#' @param data vector of numeric data
#' @param method one of "expected" and "observed", calculating the expected / observed information matrix
#'
#' @name infoMat_ast
#' @examples
#' pars <- c(mu = 0.12, sigma = 0.6, alpha = 0.7, nu1 = 3, nu2 = 5)
#' data <- rast(1000, 0.12, 0.6, 0.3, 3, 5)
#' infoMat <- infoMat_ast(pars, data, "expected")

#' @rdname infoMat
#' @export
infoMat_ast <- function(pars, data = c(), method = c("expected", "observed")) {
  method <- match.arg(method, c("expected", "observed"))
  mu <- pars["mu"]
  sigma <- pars["sigma"]
  alpha <- pars["alpha"]
  nu1 <- pars["nu1"]
  nu2 <- pars["nu2"]

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
    if (length(data) == 0) stop("data must be provided if 'observed' method is used")
    y <- data
    y1 <- y[y <= mu]
    y2 <- y[y > mu]

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

#' @rdname infoMat
#' @export
infoMat_sst <- function(pars) {
  mu <- pars["mu"]
  sigma <- pars["sigma"]
  alpha <- pars["alpha"]
  nu <- pars["nu"]

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
