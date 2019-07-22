#' @title Methods for astfit class
#'
#' @description Methods for astfit class,
#' this is not very well documented at the moment
#'
#' @param fit A AST fit object of class \code{\link{astfit}}
#' @param n the n-th moment to be calculated (needs better expression)
#' @param method one of "numerical" and "analytical", calculating the moments using numerical integration / analytical formula
#' only centralized moments can be calculated by analytical formula, while all moments can be calculated by numerical integration
#'
#' @name astfit-methods
#' @examples
#' data <- rast(1000, 0.12, 0.6, 0.3, 3, 5)
#' solver_control <- list('algorithm' = 'NLOPT_LN_BOBYQA', 'maxeval' = 1.0e5, 'xtol_rel' = 1.0e-8)
#' fit <- astfit(data, solver = 'nloptr', solver_control = solver_control)
#' summary(fit)
#' moment(fit, 1)
#' plot(fit)

#' @rdname astfit-methods
#' @export
summary.astfit <- function(fit) {
  fit$data <- NULL
  fit
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

#' @rdname astfit-methods
#' @export
moments.astfit <- function(fit, n, method = c("analytical", "numerical")) {
  pars <- fit$fitted_pars
  c(mean = mean_ast(pars = pars, method = method),
    variance = var_ast(pars = pars, method = method),
    skewness = skew_ast(pars = pars, method = method),
    kurtosis = kurt_ast(pars = pars, method = method))
}

# putting it here temporarily
#' @export
moments <- function(x, ...) {
  UseMethod("moments", x)
}
