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
  dist <- "GAT"
  pars <- rbind(fit$start_pars, fit$fixed_pars)
  res <- rbind(fit$fitted_pars, fit$standard_errors)
  colnames(pars) <- colnames(res) <- names(fit$fitted_pars)
  rownames(pars) <- c("start_pars", "fixed_pars")
  rownames(res) <- c("fitted_pars", "standard_errors")

  cat("Distribution: ", dist, "\n")
  cat("Observations: ", length(fit$data), "\n")
  cat("\nResult:\n")
  print(res)
  cat("\nLog-likelihood", fit$objective)
  cat("\n\nSolver: ", fit$solver)
  cat("\n\n")
  print(pars)
  cat("\nTime elapsed: ", fit$time_elapsed)
  cat("\nConvergence Message: ", fit$message)
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
  dist <- "GAT"
  res <- rbind(fit$fitted_pars, fit$standard_errors)
  colnames(res) <- names(fit$fitted_pars)
  rownames(res) <- c("fitted_pars", "standard_errors")

  cat("Distrifitbution: ", dist, "\n")
  cat("Observations: ", length(fit$data), "\n")
  cat("\nResult:\n")
  print(res)
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
