#' @title AST plots
#'
#' @description Plot methods for the AST fit class
#' @param fit A AST fit object of class \code{\link{astfit}}
#' @name plot-methods

#' @rdname plot-methods
#' @export
plot_density = function(fit) {
  list2env(fit$fitted.pars, envir = parent.frame())
  data = fit$data
  hist(data, breaks = 50, prob = TRUE)
  par(new = TRUE)
  # lines(density(data))
  x = seq(min(data), max(data), length.out = 1000)
  y = dast(x, mu, sigma, alpha, nu1, nu2)
  plot(x, y, axes = FALSE, xlab= "", ylab = "", col = 4, type = "l")
}

#' @rdname plot-methods
#' @export
qqplot = function(fit, method = "normal") {
  y = fit$data
  list2env(fit$fitted.pars, envir = parent.frame())
  p = past(y, mu, sigma, alpha, nu1, nu2)
  if(method == "normal") {
    x = qnorm(p, mean = mu, sd = sigma)
  } else if (method == "t") {
    x = qt(p, df = max(nu1, nu2))
  } else if (method == "ast") {
    x = qast(p, mu, sigma, alpha, nu1, nu2)
  }
  qqplot(x, y)
  abline(0, 1, col = 2, lty = 2)
}
