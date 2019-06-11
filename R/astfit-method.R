#' @title AST fit
#' the object currently have a S3 structure which may be shifted into S4 in future development
#' @description Method for fitting a AST distribution
#' @param spec An AST specification object of class \code{\link{astspec}}
#' @param solver Optimizer used for fitting, one of "nloptr", ...
#' @param solver_control Control arguments list passed to the optimizer.
#' @name astfit
#' @examples
#' data = rast(1000, 1.5, 2, 0.8, 3, 4)
#' spec = astspec(data)
#' solver_control = list("algorithm" = "NLOPT_LN_COBYLA", "maxeval" = 1.0e5, "xtol_rel" = 1.0e-8)
#' fit = astfit(spec, "nloptr", solver_control)

#' @rdname astfit
#' @export
# fit function for ast distribution
astfit = function(spec, solver, solver_control) {
  if(class(spec) != "astspec") stop("spec must be an astspec object")

  lb = c(-Inf, 0, 0, 0, 0)
  ub = c(Inf, Inf, 1, Inf, Inf)

  # will be expanded by number of solvers developed
  # should also implement a fixed parameter version for this
  if (solver == "nloptr") {
    res =  nloptr::nloptr(x0 = spec$start_pars,
                  eval_f = llast,
                  lb = lb,
                  ub = ub,
                  opts = solver_control,
                  y = data)
  }

  # this is temporary, fitted should be a list with its own elements
  sol_res = res #list(res$pars, ...)
  fitted_pars = res$solution

  structure(list(data = spec$data,
                 start_pars = spec$start_pars, fixed_pars = spec$fixed_pars, fitted_pars = fitted_pars,
                 solver = solver, solver_control = solver_control, sol_res = sol_res),
            class = "astfit")
}

# log-likelihood function of the AST distributions
# pars: parameter values
# y: data which you fit the distribution on
llast = function(pars, y) {
  mu = pars[1]; sigma = pars[2]; alpha = pars[3]; nu1 = pars[4]; nu2 = pars[5]
  T_ = length(y)
  y1 = y[y <= mu]
  y2 = y[y > mu]
  logl = -T_ * log(sigma) - 0.5 * (nu1 + 1) * sum(log(1 + ((y1 - mu)/(2 * alpha * sigma * K(nu1)))^2/nu1)) -
    0.5 * (nu2 + 1) * sum(log(1 + ((y2 - mu)/(2 * (1-alpha) * sigma * K(nu2)))^2/nu2))
  -logl
}
