#' @title Information Matrix Function of Asymmetric Student-t distribution
#'
#' @name ast-infomat
#' @aliases astInfoMat
#'
#' @description Information matrix, asymptotic covariance and correlation matrix functions of Asymmetric Student-t distribution
#'
#' @param pars a vector of parameter values for an AST distribution
#' @param data a vector of numeric data used to calculate observed information matrix
#' @param method one of "expected" and "observed", calculating the expected / observed information matrix
#'
#' @details
#' The expected information matrix is calculated by the expectation of the outer product of score functions,
#' analytical formulas are provided in \emph{Zhu and Galbraith(2010)}.
#' The observed information matrix is calculated by the expectation of negative Hessian Matrix of the log-likelihood function.
#'
#' @references
#' Zhu, D., & Galbraith, J. W. (2010). A generalized asymmetric Student-t distribution with application to financial econometrics. Journal of Econometrics, 157(2), 297-305.
#' \url{https://doi.org/10.1016/j.jeconom.2010.01.013}
#'
#' @examples
#' pars <- c(0.12, 0.6, 0.6, 3, 5)
#' data <- rast(1000, pars = pars)
#' 
#' round(astInfoMat(pars, data, "observed"), 4)
#' round(astInfoMat(pars, "expected"), 4)
#' 
#' round(astCov(pars), 4)
#' round(astCor(pars), 4)
#' 
#' @importFrom stats cov2cor

#' @rdname ast-infomat
#' @export
astInfoMat <- function(pars, data = c(), method = c("expected", "observed")) {
  method <- match.arg(method, c("expected", "observed"))
  mu <- pars[1]
  sigma <- pars[2]
  alpha <- pars[3]
  nu1 <- pars[4]
  nu2 <- pars[5]

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
    r1 <- length(y1)/length(y)
    r2 <- length(y2)/length(y)

    # outer product of score
    S_mu2 <- mean( c(1/nu1 * ((nu1+1)/(2*alpha*sigma*K(nu1)))^2 * (1/L(pars, y1) - 1/L(pars, y1)^2) ,
      1/nu2 * ((nu2+1)/(2*(1-alpha)*sigma*K(nu2)))^2 * (1/R(pars, y2) - 1/R(pars, y2)^2)) )
    S_sigma2 <- -1/sigma^2 + mean( c((nu1+1)/sigma^2 * (1 - 1/L(pars, y1))*(1 + 2/L(pars, y1)),
      (nu2+1)/sigma^2 * (1 - 1/R(pars, y2))*(1 + 2/R(pars, y2))) )
    S_alpha2 <- mean( c( (nu1+1)/alpha^2 * (1 + 1/L(pars, y1) - 2/L(pars, y1)^2),
      (nu2+1)/(1-alpha)^2 * (1 + 1/R(pars, y2) - 2/R(pars, y2)^2) ) )
    S_nu12 <- (- ( D(nu1) + (nu1 + 1)/2*Dprime(nu1) ) * mean( 1 - 1/L(pars, y1) ) +
      (nu1+1)/2*D(nu1)^2 * mean( 1/L(pars, y1) - 1/L(pars, y1)^2 )) * r1
    S_nu22 <- (- ( D(nu2) + (nu2 + 1)/2*Dprime(nu2) ) * mean( 1 - 1/R(pars, y2) ) +
      (nu2+1)/2*D(nu2)^2 * mean( 1/R(pars, y2) - 1/R(pars, y2)^2 )) * r2

    S_musigma <- mean( c( (nu1+1)/sigma * 1/L(pars, y1)^2 / nu1*2*(y1-mu)/(2*alpha*sigma*K(nu1))^2,
      (nu2+1)/sigma * 1/R(pars, y2)^2 / nu2*2*(y2-mu)/(2*(1-alpha)*sigma*K(nu2))^2 ) )
    S_mualpha <- mean( c( (nu1+1)/alpha * 1/L(pars, y1)^2 / nu1*2*(y1-mu)/(2*alpha*sigma*K(nu1))^2,
      -(nu2+1)/(1-alpha) * 1/R(pars, y2)^2 / nu2*2*(y2-mu)/(2*(1-alpha)*sigma*K(nu2))^2 ) )
    # mu nu1 may not be accurate
    S_munu1 <- (-1/2 * mean( 1/L(pars, y1)/nu1 * 2*(y1-mu)/(2*alpha*sigma*K(nu1))^2 ) +
      (nu1+1)/2*D(nu1) * mean( 1/L(pars, y1)^2/nu1 * 2*(y1-mu)/(2*alpha*sigma*K(nu1))^2 )) * r1
    S_munu2 <- (-1/2 * mean( 1/R(pars, y2) / nu2 * 2*(y2-mu)/(2*(1-alpha)*sigma*K(nu2))^2 ) +
      (nu2+1)/2*D(nu2) * mean( 1/R(pars, y2)^2 / nu2 * 2*(y2-mu)/(2*(1-alpha)*sigma*K(nu2))^2 )) * r2

    S_sigmaalpha <- mean( c( 2*(nu1+1)/(sigma*alpha) * ( 1/L(pars, y1) - 1/L(pars, y1)^2 ),
      -2*(nu2+1)/(sigma*(1-alpha)) * ( 1/R(pars, y2) - 1/R(pars,y2)^2 ) ) )
    S_sigmanu1 <- (-1/sigma * mean( 1 - 1/L(pars, y1) ) +
      (nu1+1)/sigma * D(nu1) * mean( 1/L(pars, y1) - 1/L(pars, y1)^2 )) * r1
    S_sigmanu2 <- (-1/sigma * mean( 1 - 1/R(pars, y2) ) +
      (nu2+1)/sigma * D(nu2) * mean( 1/R(pars, y2) - 1/R(pars, y2)^2 )) * r2

    S_alphanu1 <- (-1/alpha * mean(1 - 1/L(pars, y1)) +
      (nu1+1)/alpha * D(nu1) * mean(1/L(pars, y1) - 1/L(pars, y1)^2)) * r1
    # alpha nu2
    S_alphanu2 <- (1/(1-alpha) * mean(1 - 1/R(pars, y2)) -
      (nu2+1)/(1-alpha) * D(nu2) * mean(1/R(pars, y2) - 1/R(pars, y2)^2)) * r2
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

#' @rdname ast-infomat
#' @export
astCov <- function(pars) {
  solve(astInfoMat(pars))
}

#' @rdname ast-infomat
#' @export
astCor <- function(pars) {
  cov2cor(astCov(pars))
}
