#' @title Methods for astfit class
#'
#' @description Methods for astfit class
#'
#' @param fit A AST fit object of class \code{\link{astfit}}
#' @param method one of "numerical" and "analytical", calculating the moments using numerical integration / analytical formula
#' @param type one of "density" and "QQplot"
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
  fit$solver_control <- NULL
  fit$solver <- NULL
  fit$objective <- NULL
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
plot.astfit <- function(fit, type = NULL, dist = "ast", ...) {
  if (is.null(type)) {
    selection <- 1
    while (selection) {
      selection <- menu(c("Density", "qqplot"), title = "Make a plot selection (or 0 to exit)")
      if (selection == 1) {
        density_ast(fit, ...)
      } else if(selection == 2) {
        qqplot_ast(fit, dist, ...)
      }
    }
  } else {
    if (type == "density") {
      density_ast(fit, ...)
    } else if(type == "qqplot") {
      qqplot_ast(fit, dist, ...)
    }
  }
}

#' @rdname astfit-methods
#' @export
print.astfit <- function(fit) {
  fit$data <- NULL
  fit$sol_res <- NULL
  fit$solver_control <- NULL
  fit$solver <- NULL
  for (i in 1:length(fit)) {
    print(fit[i])
  }
}
