#' @title AST plots
#'
#' @description Plot methods for the AST fit class
#' @param fit A AST fit object of class \code{\link{astfit}}
#' @name plot-methods
#' @examples
#' data = rast(1000, 1.5, 2, 0.8, 3, 4)
#' spec <- astspec(data)
#' solver_control <- list('algorithm' = 'NLOPT_LN_COBYLA', 'maxeval' = 1.0e5, 'xtol_rel' = 1.0e-8)
#' fit <- astfit(spec, 'nloptr', solver_control)
#' plot(fit)

#' @rdname plot-methods
#' @export
plot.astfit <- function(fit) {
  selection <- 1
  while (selection) {
    selection <- menu(c("Density", "QQplot"), title = "Make a plot selection (or 0 to exit)")
    if (selection == 1) {
      density_ast(fit)
    } else if(selection == 2) {
      qqplot_ast(fit, dist = "normal")
    }
  }
}

density_ast <- function(fit) {
  data <- fit$data
  pars <- fit$fitted_pars
  mu <- pars["mu"]
  sigma <- pars["sigma"]
  alpha <- pars["alpha"]
  nu1 <- pars["nu1"]
  nu2 <- pars["nu2"]

  hist(data, breaks = 50, prob = TRUE)
  par(new = TRUE)
  # lines(density(data))
  x <- seq(min(data), max(data), length.out = 1000)
  y <- dast(x, mu, sigma, alpha, nu1, nu2)
  plot(x, y, axes = FALSE, xlab = "", ylab = "", col = 4, type = "l")
  abline(v = mu, col = 2)
}

qqplot_ast <- function(fit, dist = "normal", ...) {
  y <- fit$data
  pars <- fit$fitted_pars
  mu <- pars["mu"]
  sigma <- pars["sigma"]
  alpha <- pars["alpha"]
  nu1 <- pars["nu1"]
  nu2 <- pars["nu2"]

  p <- past(y, mu, sigma, alpha, nu1, nu2)
  if (dist == "normal") {
    x <- qnorm(p, mean = mean(y), sd = sigma)
    #px <- qnorm(c(0.25, 0.75), mean = mu, sd = sigma)
  } else if (dist == "t") {
    df = 0.5 * (nu1 + nu2)
    x <- qt(p, df = df)
    #px <- qast(c(0.25, 0.75), mu, sigma, 0.5, df, df)
  } else if (dist == "ast") {
    x <- qast(p, mu, sigma, alpha, nu1, nu2)
    #px <- qast(c(0.25, 0.75), mu, sigma, alpha, nu1, nu2)
  }
  qqplot(x, y, ...)
  px <- quantile(x, c(0.25, 0.75))
  py <- quantile(y, c(0.25, 0.75))
  slope <- diff(py)/diff(px)
  int <- py[1L] - slope * px[1L]
  abline(int, slope)
  points(px, py, col = 2)
}

