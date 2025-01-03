#' @importFrom graphics hist lines points
#' @importFrom stats ppoints qnorm quantile

density_gat <- function(fit, main = "Histogram/Density plot", ...) {
  data <- fit$data
  pars <- fit$fitted_pars

  hist(data, breaks = 50, prob = TRUE, main = main, ...)
  x <- seq(min(data), max(data), length.out = 1000)
  y <- dgat(x, pars = pars)
  lines(x, y, xlab = "", ylab = "", col = 4)
  # abline(v = mu, col = 2)
}

qqplot_gat <- function(fit, dist = "gat", main = "QQPlot", envelope = 0.95, ...) {
  y <- as.numeric(fit$data)
  pars <- fit$fitted_pars

  y <- y[order(y)]
  p <- ppoints(length(y))
  #p <- pgat(y, mu, sigma, alpha, nu1, nu2)
  if (dist == "normal") {
    x <- qnorm(p, mean = mean(y), sd = pars["phi"])
    px <- qnorm(c(0.25, 0.75), mean = pars["mu"], sd = pars["phi"])
  } else if (dist == "gat") {
    x <- qgat(p, pars = pars)
    px <- qgat(c(0.25, 0.75), pars = pars)
  } else {
    stop("dist must be one of normal and gat")
  }
  plot(x, y, main = main,
       xlab = paste(dist, "distribution"), ylab = "empirical distribution", ...)
  py <- quantile(y, c(0.25, 0.75))
  slope <- diff(py)/diff(px)
  int <- py[1L] - slope * px[1L]
  abline(int, slope, col = 4)
  points(px, py, col = 2)

  if ( envelope != FALSE) {
    zz <- qnorm(1 - (1 - envelope)/2)
    #eval(parse(text=paste("	SE <- (b/d.function(z,", distributionParameter,",...))* sqrt(P * (1 - P)/n)")))
    SE <- slope/dgat(x, pars = pars) * sqrt(p * (1 - p)/length(y))
    #		SE <- (b/d.function(z, ...)) * sqrt(P * (1 - P)/n)
    fit.value <- int + slope * x
    upper <- fit.value + zz * SE
    lower <- fit.value - zz * SE
    lines(x, upper, lty = 2, col = 4)
    lines(x, lower, lty = 2, col = 4)
  }
}
