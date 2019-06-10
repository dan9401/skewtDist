#' @title AST Specification
#' the object currently have a S3 structure which will be shifted into S4 in future development
#' @description Constructor of a AST distribution object before fitting the parameters
#' @param start.pars List of starting parameters for the optimization algorithm
#' @param fixed.pars List of parameters to be kept fixed during the optimization routine.
#' @name astspec

#' @rdname astspec
#' @export
astspec = function(start.pars = list(), fixed.pars = list()) {
  # if(!is.list(start.pars)) stop("start.pars must be a named list")
  # if(!is.list(fixed.pars)) stop("fixed.pars must be a named list")
  # do.call(check_bound, start.pars)
  # do.call(check_bound, fixed.pars)
  structure(list(start.pars = start.pars, fixed.pars = fixed.pars),
            class = "astspec")
}

# check boundary
check_bound = function(mu, sigma, alpha, nu1, nu2) {
  if(!is.numeric(mu)) stop("mu must be numeric")
  if(!is.numeric(sigma)) stop("mu must be numeric")
  if(!is.numeric(alpha)) stop("mu must be numeric")
  if(!is.numeric(nu1)) stop("mu must be numeric")
  if(!is.numeric(nu2)) stop("mu must be numeric")
  if(sigma <= 0) stop("sigma must be greater than 0")
  if(nu1 <= 0) stop("nu1 must be greater than 0")
  if(nu2 <= 0) stop("nu2 must be greater than 0")
  if(alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
}
