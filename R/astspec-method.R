#' @title AST Specification
#' the object currently have a S3 structure which will be shifted into S4 in future development
#' @description Constructor of a AST distribution object before fitting the parameters
#' @param start.pars List of starting parameters for the optimization algorithm
#' @param fixed.pars List of parameters to be kept fixed during the optimization routine.
#' @name astspec

#' @rdname astspec
#' @export
astspec = function(start.pars = list(), fixed.pars = list()) {
  if(!is.list(start.pars)) stop("start.pars must be a named list")
  if(!is.list(fixed.pars)) stop("fixed.pars must be a named list")
  do.call(check_bound, start.pars)
  do.call(check_bound, fixed.pars)
  structure(list(start.pars = start.pars, fixed.pars = fixed.pars),
            class = "astspec")
}
