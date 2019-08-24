#' @title Methods for gat class
#'
#' @description Methods for gat S3 class
#'
#' @param fit A GAT fit object of class \code{\link{gat}}
#' @param method one of "numerical" and "analytical", calculating the moments using numerical integration / analytical formula
#' @param type one of "density" and "QQplot"
#'
#' @details should also add the empirical moments
#'
#' @name gat-methods
#' @aliases summary.gat
#' @aliases moments.gat
#' @aliases plots.gat
#' @aliases fitted.gat
#' @aliases se.gat
#' @aliases objective.gat
#'
#' @examples
#' pars <- c(0.12, 0.6, 0.6, 6, 5)
#' data <- rgat(1000, pars = pars)
#' solver_control <- list(eval.max = 10^3, iter.max = 10^3)
#' fit <- gatMLE(data, solver = 'nlminb', solver_control = solver_control)
#' summary(fit)
#' moments(fit)
#' fitted(fit)
#' se(fit)
#' objective(fit)
#' plot(fit)

#' @rdname gat-methods
#' @export
summary.gat <- function(fit) {
  fit$data <- NULL
  fit$solver_control <- NULL
  fit$solver <- NULL
  fit$objective <- NULL
  fit
}

#' @rdname gat-methods
#' @export
moments.gat <- function(fit, method = c("analytical", "numerical")) {
  pars <- fit$fitted_pars
  gatMoments(pars = pars, method)
}

#' @rdname gat-methods
#' @export
print.gat <- function(fit) {
  fit$data <- NULL
  fit$sol_res <- NULL
  fit$solver_control <- NULL
  fit$solver <- NULL
  for (i in 1:length(fit)) {
    print(fit[i])
  }
}

#' @rdname gat-methods
#' @export
fitted.gat <- function(fit) {
  fit$fitted_pars
}

#' @rdname gat-methods
#' @export
se.gat <- function(fit) {
  fit$standard_errors
}

#' @rdname gat-methods
#' @export
objective.gat <- function(fit) {
  fit$objective
}

#' @rdname gat-methods
#' @export
plot.gat <- function(fit, type = NULL, dist = "gat", envelope = 0.95, ...) {
  if (is.null(type)) {
    selection <- 1
    while (selection) {
      selection <- menu(c("Density", "qqplot"), title = "Make a plot selection (or 0 to exit)")
      if (selection == 1) {
        density_gat(fit, ...)
      } else if(selection == 2) {
        qqplot_gat(fit, dist, envelope = 0.95, ...)
      }
    }
  } else {
    if (type == "density") {
      density_gat(fit, ...)
    } else if(type == "qqplot") {
      qqplot_gat(fit, dist = dist, envelope = envelope, ...)
    }
  }
}
