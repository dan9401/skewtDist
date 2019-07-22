# Plot methods for the AST and GAT distribuitons, the aesthetics would require further adjustment, either adjust the graphic parameters, or a new graphics engine
# we may also want to keep separate files for both distributions, doesn't seem necessary at the time
# and plot method for ast & gat class without data, just for exploration uses
# authorized domain, here or in the moment file

density_ast <- function(fit, ...) {
  data <- fit$data
  pars <- fit$fitted_pars
  mu <- pars["mu"]
  sigma <- pars["sigma"]
  alpha <- pars["alpha"]
  if (fit$symmetric == TRUE) {
    nu1 <- nu2 <- pars["nu"]
  } else {
    nu1 <- pars["nu1"]
    nu2 <- pars["nu2"]
  }

  hist(data, breaks = 50, prob = TRUE)
  par(new = TRUE)
  # lines(density(data))
  x <- seq(min(data), max(data), length.out = 1000)
  y <- dast(x, mu, sigma, alpha, nu1, nu2)
  plot(x, y, axes = FALSE, xlab = "", ylab = "", col = 4, type = "l", ...)
  abline(v = mu, col = 2)
}

qqplot_ast <- function(fit, dist = "ast", ...) {
  y <- fit$data
  pars <- fit$fitted_pars
  mu <- pars["mu"]
  sigma <- pars["sigma"]
  alpha <- pars["alpha"]
  if (fit$symmetric == TRUE) {
    nu1 <- nu2 <- pars["nu"]
  } else {
    nu1 <- pars["nu1"]
    nu2 <- pars["nu2"]
  }

  y <- y[order(y)]
  p <- ppoints(length(y))
  #p <- past(y, mu, sigma, alpha, nu1, nu2)
  if (dist == "normal") {
    x <- qnorm(p, mean = mean(y), sd = sigma)
    px <- qnorm(c(0.25, 0.75), mean = mu, sd = sigma)
  } else if (dist == "ast") {
    x <- qast(p, mu, sigma, alpha, nu1, nu2)
    px <- qast(c(0.25, 0.75), mu, sigma, alpha, nu1, nu2)
  }
  plot(x, y, ...)
  py <- quantile(y, c(0.25, 0.75))
  slope <- diff(py)/diff(px)
  int <- py[1L] - slope * px[1L]
  abline(int, slope, col = 4)
  points(px, py, col = 2)
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

  start_pars <- c(mu = 0, sigma = 1, alpha = 0.5, nu1 = 2, nu2 = 2)
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
