#' @title AST fit
#' the object currently have a S3 structure which will be shifted into S4 in future development
#' @description Method for fitting a AST distribution
#' @param spec An AST specification object of class \code{\link{astspec}}
#' @param data A univariate data object to fit on
#' @param solver Optimization algorithm used, one of ...
#' @param solver.control Control arguments passed to the optimization algorithm ...
#' @param fit.control Control arguments passed to the fitting routine...
#' @name astfit
#' @examples
#' spec = astspec(c(0, 1, 0.5, 1, 1))
#' data = rast(1000, 1.5, 1.2, 0.8, 3, 4)
#' fit = astfit(spec, data, "nloptr", list("algorithm" = "NLOPT_LN_COBYLA", "maxeval" = 1.0e5, "xtol_rel" = 1.0e-8))

#' @rdname astfit
#' @export
astfit = function(spec, data, solver, solver.control) {
  if(!is.numeric(data)) stop("data must be numeric")
  if(class(spec) != "astspec") stop("spec must be an astspec object")

  lb = c(-Inf, 0, 0, 0, 0)
  ub = c(Inf, Inf, 1, Inf, Inf)

  # The line below will be substituted by the optimization algorithm in the actual implementation
  # maybe expanded by number of solvers developed
  if (solver == "nloptr") {
    res =  nloptr::nloptr(x0 = spec$start.pars,
                  eval_f = llast,
                  lb = lb,
                  ub = ub,
                  opts = solver.control,
                  y = data)
  }

  fitted = res #list(res$pars, ...)
  structure(list(data = data, start.pars = spec$start.pars, fitted = fitted),
            class = "astfit")
}

# log-likelihood function of the AST distributions
# pars: parameter values
# y: data which you fit the distribution on
llast = function(pars, y) {
  #list2env(pars, envir = parent.frame())
  mu = pars[1]; sigma = pars[2]; alpha = pars[3]; nu1 = pars[4]; nu2 = pars[5];
  T_ = length(y)
  y1 = y[y <= mu]
  y2 = y[y > mu]
  logl = -T_ * log(sigma) - 0.5 * (nu1 + 1) * sum(log(1 + ((y1 - mu)/(2 * alpha * sigma * K(nu1)))^2/nu1)) -
    0.5 * (nu2 + 1) * sum(log(1 + ((y2 - mu)/(2 * (1-alpha) * sigma * K(nu2)))^2/nu2))
  -logl
}
