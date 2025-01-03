#' @title Methods for ast class
#'
#' @description Methods for ast S3 class
#'
#' @param object A AST fit object of class \code{\link{ast}}
#' @param x A AST fit object of class \code{\link{ast}}
#' @param method one of "numerical" and "analytical", calculating the moments using numerical integration / analytical formula
#' @param type one of "density" or "qqplot"
#' @param dist one of "norm" or "ast", the theoretical distribution used in QQplots
#' @param envelope the confidence level used to construct the envelope
#' @param ... additional arguments for the \code{hist} or \code{plot} function from \code{graphics}
#'
#' @details should also add the empirical moments
#'
#' @name ast-methods
#' @aliases summary.ast
#' @aliases moments.ast
#' @aliases plots.ast
#'
#' @examples
#' pars <- c(0.12, 0.6, 0.6, 6, 5)
#' data <- rast(1000, pars = pars)
#' solver_control <- list(eval.max = 10^3, iter.max = 10^3)
#' fit <- astMLE(data, solver = 'nlminb', solver_control = solver_control)
#' summary(fit)
#' moments(fit)
#' plot(fit, 1)
#' 
#' @importFrom utils menu

#' @rdname ast-methods
#' @export
summary.ast <- function(object, ...) {
  fit <- object
  dist <- ifelse(fit$symmetric == TRUE, "SST", "AST")
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
  cat("\n")
}

#' @rdname ast-methods
#' @export
moments.ast <- function(x, method = c("analytical", "numerical"), ...) {
  fit <- x
  pars <- fit$fitted_pars
  astMoments(pars = pars, method)
}

#' @rdname ast-methods
#' @export
print.ast <- function(x, ...) {
  fit <- x
  dist <- ifelse(fit$symmetric == TRUE, "SST", "AST")
  res <- rbind(fit$fitted_pars, fit$standard_errors)
  colnames(res) <- names(fit$fitted_pars)
  rownames(res) <- c("fitted_pars", "standard_errors")

  cat("Distrifitbution: ", dist, "\n")
  cat("Observations: ", length(fit$data), "\n")
  cat("\nResult:\n")
  print(res)
}

#' @rdname ast-methods
#' @export
plot.ast <- function(x, type = NULL, dist = "ast", envelope = 0.95, ...) {
  fit <- x
  if (is.null(type)) {
    selection <- 1
    while (selection) {
      selection <- menu(c("Density", "QQplot"), title = "Make a plot selection (or 0 to exit)")
      if (selection == 1) {
        density_ast(fit, ...)
      } else if(selection == 2) {
        qqplot_ast(fit, dist, envelope = 0.95, ...)
      }
    }
  } else {
    if (type == "density") {
      density_ast(fit, ...)
    } else if(type == "qqplot") {
      qqplot_ast(fit, dist = dist, envelope = envelope, ...)
    }
  }
}
