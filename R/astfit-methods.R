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

#' @rdname astfit-methods
#' @export
moment.astfit <- function(fit, n, method = c("numerical", "analytical")) {
  pars <- fit$fitted_pars
  mu <- pars["mu"]
  sigma <- pars["sigma"]
  alpha <- pars["alpha"]
  nu1 <- pars["nu1"]
  nu2 <- pars["nu2"]
  method <- match.arg(method)
  if (method == "analytical") {
    if (mu != 0) {
      stop("Analytical formula of moments cannot calculate with location parameters other than 0.")
    }
    # return value - analytical
    moment_ast_analytical(n, mu, sigma, alpha, nu1, nu2)
  } else {
    # return value - numerical
    moment_ast_numerical(n, mu, sigma, alpha, nu1, nu2)
  }
}

# putting it here temporarily
#' @export
moment <- function(x, ...) {
  UseMethod("moment", x)
}
