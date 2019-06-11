#' @title AST Specification
#' the object currently have a S3 structure which may be shifted into S4 in future development
#' @description Method for creating an AST distribution object prior to fitting.
#' @param data A univariate data object, can be ... for the AST distribution to fit on.
#' @param start_pars Numeric vector of starting parameters for the optimization algorithm.
#' @param fixed_pars Numeric vector of parameters to be kept fixed during the optimization routine.
#' @name astspec
#' @examples
#' data = rast(1000, 1.5, 2, 0.8, 3, 4)
#' spec = astspec(data)

#' @rdname astspec
#' @export
astspec = function(data, start_pars = c(0, 1, 0.5, 1, 1), fixed_pars = c()) {
  if (!is.numeric(data)) stop("data must be numeric")
  if (length(start_pars) != 5) stop("start_pars must be a numeric of length 5")
  # if (!is.list(fixed_pars)) stop("fixed_pars must be a named list")

  check_bound(start_pars)
  # if (length(fixed_pars) != 0) do.call(check_bound, as.list(fixed_pars))
  structure(list(data = data, start_pars = start_pars, fixed_pars = fixed_pars),
            class = "astspec")
}

# check boundary
check_bound = function(pars) {
  mu = pars[1]; sigma = pars[2]; alpha = pars[3]; nu1 = pars[4]; nu2 = pars[5]
  if(!is.numeric(mu)) stop("mu must be numeric")
  if(!is.numeric(sigma)) stop("sigma must be numeric")
  if(!is.numeric(alpha)) stop("alpha must be numeric")
  if(!is.numeric(nu1)) stop("nu1 must be numeric")
  if(!is.numeric(nu2)) stop("nu2 must be numeric")
  if(sigma <= 0) stop("sigma must be greater than 0")
  if(nu1 <= 0) stop("nu1 must be greater than 0")
  if(nu2 <= 0) stop("nu2 must be greater than 0")
  if(alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
}
