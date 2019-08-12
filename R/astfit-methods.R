#' @title Methods for astFit class
#'
#' @description Methods for astFit S3 class
#'
#' @param fit A AST fit object of class \code{\link{astFit}}
#' @param method one of "numerical" and "analytical", calculating the moments using numerical integration / analytical formula
#' @param type one of "density" and "QQplot"
#'
#' @details should also add the empirical moments
#'
#' @name astFit-methods
#' @aliases summary.astFit
#' @aliases moments.astFit
#' @aliases plots.astFit
#' @aliases fitted.astFit
#' @aliases se.astFit
#' @aliases objective.astFit
#'
#' @examples
#' pars <- c(0.12, 0.6, 0.6, 6, 5)
#' data <- rast(1000, pars = pars)
#' solver_control <- list(eval.max = 10^3, iter.max = 10^3)
#' fit <- astFit(data, solver = 'nlminb', solver_control = solver_control)
#' summary(fit)
#' moments(fit)
#' fitted(fit)
#' se(fit)
#' objective(fit)
#' plot(fit)

#' @rdname astFit-methods
#' @export
summary.astFit <- function(fit) {
  fit$data <- NULL
  fit$solver_control <- NULL
  fit$solver <- NULL
  fit$objective <- NULL
  fit
}

#' @rdname astFit-methods
#' @export
moments.astFit <- function(fit, method = c("analytical", "numerical")) {
  pars <- fit$fitted_pars
  astMoments(pars = pars, method)
}

#' @rdname astFit-methods
#' @export
print.astFit <- function(fit) {
  fit$data <- NULL
  fit$sol_res <- NULL
  fit$solver_control <- NULL
  fit$solver <- NULL
  for (i in 1:length(fit)) {
    print(fit[i])
  }
}

#' @rdname astFit-methods
#' @export
fitted.astFit <- function(fit) {
  fit$fitted_pars
}

#' @rdname astFit-methods
#' @export
se.astFit <- function(fit) {
  fit$standard_errors
}

#' @rdname astFit-methods
#' @export
objective.astFit <- function(fit) {
  fit$objective
}

#' @rdname astFit-methods
#' @export
plot.astFit <- function(fit, type = NULL, dist = "ast", ...) {
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
