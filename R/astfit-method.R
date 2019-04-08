#' @title AST fit
#' the object currently have a S3 structure which will be shifted into S4 in future development
#' @description Method for fitting a AST distribution
#' @param spec An AST specification object of class \code{\link{astspec}}
#' @param data A univariate data object to fit on
#' @param solver Optimization algorithm used, one of ...
#' @param solver.control Control arguments passed to the optimization algorithm ...
#' @param fit.control Control arguments passed to the fitting routine...
#' @name astfit

#' @rdname astfit
#' @export
astfit = function(spec, data, solver, solver.control, fit.control) {
  if(!is.numeric(data)) stop("data must be numeric")
  if(class(spec) != "astspec") stop("spec must be an astspec object")

  # The line below will be substituted by the optimization algorithm in the actual implementation
  fitted.pars = nlm(ll, data, spec$start.pars)
  structure(list(data = data, start.pars = spec$start.pars, fixed.pars = fixed.pars, fitted.pars = fitted.pars),
            class = "astfitc")
}

# log-likelihood function of the AST distributions
ll = function(y, start.pars) {
  list2env(start.pars, envir = parent.frame())
  T_ = length(y)
  y1 = y[y <= mu]
  y2 = y[y > mu]
  logl = -T_ * log(sigma) - 0.5 * (nu1 + 1) * sum(log(1 + ((y1 - mu)/(2 * alpha * sigma * K(nu1)))^2/nu1)) -
    0.5 * (nu2 + 1) * sum(log(1 + ((y2 - mu)/(2 * (1-alpha) * sigma * K(nu2)))^2/nu2))
  logl
}
