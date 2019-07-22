#' @title Methods for astfit class
#'
#' @description Methods for astfit class
#'
#' @param fit A AST fit object of class \code{\link{astfit}}
#' @param method one of "numerical" and "analytical", calculating the moments using numerical integration / analytical formula
#' @param selection one of 1 & 2, for plotting the density or QQPlot
#'
#' @details should also add the empirical moments
#'
#' @name astfit-methods
#' @aliases summary.astfit
#' @aliases moments.astfit
#' @aliases plots.astfit
#'
#'
#' @examples
#' pars <- c(0.12, 0.6, 0.4, 3, 5)
#' data <- rast(1000, pars = pars)
#' solver_control <- list(eval.max = 10^3, iter.max = 10^3)
#' fit <- astfit(data, solver = 'nlminb', solver_control = solver_control)
#' summary(fit)
#' moments(fit)
#' plot(fit)

#' @rdname astfit-methods
#' @export
summary.astfit <- function(fit) {
  fit$data <- NULL
  fit
}

#' @rdname astfit-methods
#' @export
moments.astfit <- function(fit, method = c("analytical", "numerical")) {
  pars <- fit$fitted_pars
  moments_ast(pars = pars, method)
}

#' @rdname astfit-methods
#' @export
plot.astfit <- function(fit, selection = NULL, ...) {
  if (is.null(selection)) {
    selection <- 1
    while (selection) {
      selection <- menu(c("Density", "QQplot"), title = "Make a plot selection (or 0 to exit)")
      if (selection == 1) {
        density_ast(fit, ...)
      } else if(selection == 2) {
        qqplot_ast(fit, dist = "ast", ...)
      }
    }
  } else {
    if (selection == 1) {
      density_ast(fit, ...)
    } else if(selection == 2) {
      qqplot_ast(fit, dist = "ast", ...)
    }
  }

}
