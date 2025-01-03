#' @importFrom graphics hist lines abline points
#' @importFrom stats ppoints qnorm quantile

density_ast <- function(fit, main = "Histogram/Density plot", ...) {
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

  hist(data, breaks = 50, prob = TRUE, main = main, ...)
  x <- seq(min(data), max(data), length.out = 1000)
  y <- dast(x, mu, sigma, alpha, nu1, nu2)
  lines(x, y, xlab = "", ylab = "", col = 4)
  abline(v = mu, col = 2)
}

qqplot_ast <- function(fit, dist = "ast", main = "QQPlot", envelope = 0.95, ...) {
  y <- as.numeric(fit$data)
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
    x <- qnorm(p, mean = 0, sd = 1)
    px <- qnorm(c(0.25, 0.75), mean = 0, sd = 1)
  } else if (dist == "ast") {
    x <- qast(p, mu, sigma, alpha, nu1, nu2)
    px <- qast(c(0.25, 0.75), mu, sigma, alpha, nu1, nu2)
  } else {
    stop("dist must be one of normal and ast")
  }

  plot(x, y, main = main,
       pch = 20, ...)
  py <- quantile(y, c(0.25, 0.75))
  slope <- diff(py)/diff(px)
  int <- py[1L] - slope * px[1L]
  abline(int, slope, col = 4)
  points(px, py, col = 2)

  if ( envelope != FALSE) {
    zz <- qnorm(1 - (1 - envelope)/2)

    #eval(parse(text=paste("	SE <- (b/d.function(z,", distributionParameter,",...))* sqrt(P * (1 - P)/n)")))

    SE <- (slope/dast(x, pars = pars)) * sqrt(p * (1 - p)/length(y))
    #		SE <- (b/d.function(z, ...)) * sqrt(P * (1 - P)/n)
    fit.value <- int + slope * x
    upper <- fit.value + zz * SE
    lower <- fit.value - zz * SE
    lines(x, upper, lty = 2, col = 4)
    lines(x, lower, lty = 2, col = 4)
  }
}
