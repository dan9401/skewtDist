#' @title AST plots
#'
#' @description Plot methods for the AST fit class
#' @param fit A AST fit object of class \code{\link{astfit}}
#' @name plot-methods
#' @examples
#' data = rast(1000, 1.5, 2, 0.8, 3, 4)
#' spec = astspec(data)
#' solver_control = list("algorithm" = "NLOPT_LN_COBYLA", "maxeval" = 1.0e5, "xtol_rel" = 1.0e-8)
#' fit = astfit(spec, "nloptr", solver_control)
#'
#' density_ast(fit)
#' qqplot_ast(fit, dist = "normal")
#' qqplot_ast(fit, dist = "t")
#' qqplot_ast(fit, dist = "ast")

#' @rdname plot-methods
#' @export
density_ast = function(fit) {
  data = fit$data
  pars = fit$fitted_pars
  mu = pars[1]; sigma = pars[2]; alpha = pars[3]; nu1 = pars[4]; nu2 = pars[5]

  hist(data, breaks = 50, prob = TRUE)
  par(new = TRUE)
  # lines(density(data))
  x = seq(min(data), max(data), length.out = 1000)
  y = dast(x, mu, sigma, alpha, nu1, nu2)
  plot(x, y, axes = FALSE, xlab= "", ylab = "", col = 4, type = "l")
}

#' @rdname plot-methods
#' @export
qqplot_ast = function(fit, dist = "normal", ...) {
  y = fit$data
  pars = fit$fitted_pars
  mu = pars[1]; sigma = pars[2]; alpha = pars[3]; nu1 = pars[4]; nu2 = pars[5]

  p = past(y, mu, sigma, alpha, nu1, nu2)
  if(dist == "normal") {
    x = qnorm(p, mean = mu, sd = sigma)
    px = qnorm(c(0.25, 0.75))
  } else if (dist == "t") {
    x = qt(p, df = max(nu1, nu2))
  } else if (dist == "ast") {
    x = qast(p, mu, sigma, alpha, nu1, nu2)
  }
  qqplot(x, y, ...)

  py = quantile(y, c(0.25, 0.75))
  slope = diff(y) / diff(x)
  int = y[1L] - slope * x[1L]
  abline(int, slope)

}

#' #' @rdname plot-methods
#' #' @export
#' plot_density.gat = function(fit) {
#'   list2env(fit$fitted_pars, envir = parent.frame())
#'   data = fit$data
#'   hist(data, breaks = 50, prob = TRUE)
#'   par(new = TRUE)
#'   # lines(density(data))
#'   x = seq(min(data), max(data), length.out = 1000)
#'   y = dast(x, mu, sigma, alpha, nu1, nu2)
#'   plot(x, y, axes = FALSE, xlab= "", ylab = "", col = 4, type = "l")
#' }
#'
#' #' @rdname plot-methods
#' #' @export
#' qqplot.gatfit = function(fit, method = "normal") {
#'   y = fit$data
#'   list2env(fit$fitted_pars, envir = parent.frame())
#'   p = past(y, mu, sigma, alpha, nu1, nu2)
#'   if(method == "normal") {
#'     x = qnorm(p, mean = mu, sd = sigma)
#'   } else if (method == "t") {
#'     x = qt(p, df = max(nu1, nu2))
#'   } else if (method == "ast") {
#'     x = qast(p, mu, sigma, alpha, nu1, nu2)
#'   }
#'   qqplot(x, y)
#'   abline(0, 1, col = 2, lty = 2)
#' }
