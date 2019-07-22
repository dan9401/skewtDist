# Plot methods for the AST and GAT distribuitons, the aesthetics would require further adjustment,
# either adjust the graphic parameters, or a new graphics engine
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
  x <- seq(min(data), max(data), length.out = 1000)
  y <- dast(x, mu, sigma, alpha, nu1, nu2)
  lines(x, y, xlab = "", ylab = "", col = 4, ...)
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
