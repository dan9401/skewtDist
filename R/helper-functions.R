# may keep separate files for gat & ast distribuitons

# these 4 are quality of life functions, to keep the formulas simple, the list may be extended
# as the current hessian is still very complicated for ast, not mentioning the unimplemented gat
K <- function(nu) {
  # there is a precision error here, numerical expert needed may need a c++ version
  # nu = 1e10 and greater start to be varying
  if (nu < 1000) {
    exp(lgamma(0.5 * (nu + 1)) - log(sqrt(pi * nu)) - lgamma(0.5 * nu))
  } else {
    # seems to be correct, need to check
    dnorm(0)
  }
}

D <-  function(nu) {
  D <- digamma((nu + 1) / 2) - digamma(nu / 2)
  D
}

Dprime <- function(nu) {
  Dp <- trigamma((nu+1)/2) - trigamma(nu/2)
}

L <- function(pars, y) {
  mu <- pars[1]
  sigma <- pars[2]
  alpha <- pars[3]
  nu1 <- pars[4]
  L <- 1 + 1 / nu1 * ( (y - mu) / (2 * alpha * sigma * K(nu1)) )^2
  L
}

R <- function(pars, y) {
  mu <- pars[1]
  sigma <- pars[2]
  alpha <- pars[3]
  nu2 <- ifelse(length(pars) == 4, pars[4], pars[5])
  R <- 1 + 1 / nu2 * ( (y - mu) / (2 * (1 - alpha) * sigma * K(nu2)) )^2
  R
}

# the following 2 are required for surface plot, for exploraton reason
parVec <- function(x, xName) {
  eps <- 1.0e-8
  if (xName == "mu") {
    seq(x - 0.5, x + 0.5, 0.1)
  } else if (xName == "alpha") {
    seq(0.1, 0.9, 0.1)
  } else {
    seq(max(0.1, x - 0.5), x + 0.5, 0.1)
  }
}

# this may be discarded or put into another form, not sure at the time
# but most likely needed
check_bound <- function(pars) {
  mu <- pars["mu"]
  sigma <- pars["sigma"]
  alpha <- pars["alpha"]
  nu1 <- pars["nu1"]
  nu2 <- pars["nu2"]

  if (!is.numeric(mu))
    stop("mu must be numeric")
  if (!is.numeric(sigma))
    stop("sigma must be numeric")
  if (!is.numeric(alpha))
    stop("alpha must be numeric")
  if (!is.numeric(nu1))
    stop("nu1 must be numeric")
  if (!is.numeric(nu2))
    stop("nu2 must be numeric")
  if (sigma <= 0)
    stop("sigma must be greater than 0")
  if (nu1 <= 0)
    stop("nu1 must be greater than 0")
  if (nu2 <= 0)
    stop("nu2 must be greater than 0")
  if (alpha <= 0 || alpha >= 1)
    stop("alpha must be between 0 and 1")
}

# used only for surfaceplot in plot-methods
obj_surface <- function(pars, data, start_pars, fixed_pars, solver, solver_control, xName, yName, symmetric = FALSE) {
  fixed_pars[xName] <- pars[1]
  fixed_pars[yName] <- pars[2]
  #astfit_local(data, start_pars, fixed_pars, solver, solver_control, symmetric)
  return(astfit_local(data, start_pars, fixed_pars, solver, solver_control, symmetric)$objective)
}

safeIntegrate <- function(f, lower, upper, ..., subdivisions = 1000L,
                          rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
                          stop.on.error = TRUE, keep.xy = FALSE, aux = NULL) {
  res <- integrate(f, lower, upper, ..., subdivisions = subdivisions,
                   rel.tol = rel.tol, abs.tol = abs.tol,
                   stop.on.error = stop.on.error, keep.xy = keep.xy, aux = aux)
  # the safeIntegrate in HyperbolicDist has different rules
  # don't seem to be necessary
  if (lower == upper) {
    res$value <- 0
    res$abs.error <- 0
  }
  res
}

newton_raphson <- function(f, x0 = 0.5, maxiter = 1e2, tol = 1e-8) {
  h <- 1e-8
  i <- 1
  x1 <- x0
  # p <- numeric(maxiter)
  while(i <= maxiter) {
    fprime <- (f(x0 + h) - f(x0)) / h
    x1 <- (x0 - (f(x0) / fprime))
    # p = x1
    i <- i + 1
    if (abs(x1 - x0) < tol) {
      break
    }
    x0 <- x1
  }
  x1
}


# putting it here temporarily
#' @export
moments <- function(x, ...) {
  UseMethod("moments", x)
}

# has no documentation developed yet
#' @export
surfacePlot <- function(n, pars, plotPars, ...) {
  mu <- pars[1]
  sigma <- pars[2]
  alpha <- pars[3]
  nu1 <- pars[4]
  nu2 <- pars[5]
  data <- rast(n, mu, sigma, alpha, nu1, nu2)

  xName <- plotPars[1]
  yName <- plotPars[2]
  xVec <- parVec(pars[xName], xName)
  yVec <- parVec(pars[yName], yName)
  xLen <- length(xVec)
  yLen <- length(yVec)
  xMat <- matrix(rep(xVec, yLen), xLen, yLen)
  yMat <- matrix(rep(yVec, xLen), xLen, yLen, byrow = TRUE)
  parGrid <- array(c(xMat, yMat), c(xLen, yLen, 2))

  start_pars <- c()
  fixed_pars <- c()
  solver <- "Rsolnp"
  solver_control <- list(trace = 0)
  valGrid <- apply(parGrid, 1:2, obj_surface, data, start_pars, fixed_pars, solver, solver_control, xName, yName)
  rownames(valGrid) <- xVec
  colnames(valGrid) <- yVec
  persp(xVec, yVec, valGrid, xlab = xName, ylab = yName, ...)
  return(list(xVec, yVec, valGrid))
  # p <- plot_ly(z = ~valueGrid) %>% add_surface()
  # chart_link = api_create(p, filename = paste("grid", plotPars[1], plotPars[2], sep = ""))
  # chart_link
}

report <- function(fit) {
  report <- c(round(c(fit$fitted_pars, fit$objective, fit$time_elapsed) , 8), fit$message)
  names(report)[c(6, 7, 8)] <- c("objective", "time", "message")
  report
}

# has no documentation developed yet
#' @export
my_report <- function(data, solver, solver_control, plots = "none", dist = "ast") {
  fitList <- lapply(data, astfit, solver = solver, solver_control = solver_control)

  if (plots == "both") {
    par(mfrow = c(4, 5))
    lapply(fitList, plot, type = "density")
    par(mfrow = c(4, 5))
    lapply(fitList, plot, type = "qqplot", dist = dist)
  #}
  } else if (plots == "density") {
    par(mfrow = c(4, 5))
    lapply(fitList, plot, type = "density")
  } else if (plots == "qqplot") {
    par(mfrow = c(4, 5))
    lapply(fitList, plot, type = "qqplot", dist = dist)
  }
  par(mfrow = c(1, 1))

  res = t(sapply(fitList, report))
  res
}
