#' @title AST fit
#' the object currently have a S3 structure which may be shifted into S4 in future development
#' @description Method for fitting a AST distribution
#' @param spec An AST specification object of class \code{\link{astspec}}
#' @param solver Optimizer used for fitting, one of 'nloptr', ...
#' @param solver_control Control arguments list passed to the optimizer.
#' @name astfit
#' @examples
#' data <- rast(1000, 1.5, 2, 0.8, 3, 4)
#' spec <- astspec(data)
#' solver_control <- list('algorithm' = 'NLOPT_LN_COBYLA', 'maxeval' = 1.0e5, 'xtol_rel' = 1.0e-8)
#' fit <- astfit(spec, 'nloptr', solver_control)

#' @rdname astfit
#' @export
# fit function for ast distribution
astfit <- function(spec, solver, solver_control) {
    if (class(spec) != "astspec")
        stop("spec must be an astspec object")

    data <- spec$data
    ipars <- spec$ipars
    est_idx <- which(is.na(ipars$fixed_pars))
    start_pars <- ipars$start_pars[est_idx]
    names(start_pars) <- rownames(ipars)[est_idx]
    fixed_pars <- ipars$start_pars[-est_idx]
    names(fixed_pars) <- rownames(ipars)[-est_idx]

    if ("nu1" %in% names(fixed_pars)) { nu1vec <- ipars["nu1", "fixed_pars"] } else { nu1vec <- seq(2, 20, by = 4) }
    if ("nu2" %in% names(fixed_pars)) { nu2vec <- ipars["nu2", "fixed_pars"] } else { nu2vec <- seq(2, 20, by = 4) }
    grid <- array(c(rep(nu1vec, length(nu2vec)), rep(nu2vec, length(nu1vec))), c(length(nu1vec), length(nu2vec), 2))
    if (dim(grid)[1] == dim(grid)[2]) { grid[,,2] = t(grid[,,2]) }
    valueGrid <- apply(grid, 1:2, objective_value, spec = spec, solver = solver, solver_control = solver_control)
    idx <- as.numeric(which(valueGrid == min(valueGrid), arr.ind = T))
    idx <- c(grid[idx[1], idx[2], 1], grid[idx[1], idx[2], 2])

    if (!("nu1" %in% names(fixed_pars))) { nu1vec <- seq(idx[1]-2, idx[1]+2, by=1) }
    if (!("nu2" %in% names(fixed_pars))) { nu2vec <- seq(idx[2]-2, idx[2]+2, by=1) }
    grid <- array(c(rep(nu1vec, length(nu2vec)), rep(nu2vec, length(nu1vec))), c(length(nu1vec), length(nu2vec), 2))
    if (dim(grid)[1] == dim(grid)[2]) { grid[,,2] = t(grid[,,2]) }
    valueGrid <- apply(grid, 1:2, objective_value, spec = spec, solver = solver, solver_control = solver_control)
    idx <- as.numeric(which(valueGrid == min(valueGrid), arr.ind = T))
    idx <- c(grid[idx[1], idx[2], 1], grid[idx[1], idx[2], 2])

    spec <- astspec(data, start_pars = start_pars, fixed_pars = c(fixed_pars[!(names(fixed_pars) %in% c("nu1", "nu2"))], "nu1" = idx[1], "nu2" = idx[2]))
    fit <- astfit_local(spec, solver, solver_control)

    start_pars <- c(fit$sol_res$solution, c(idx[1], idx[2]))
    names(start_pars) <- c("alpha", "mu", "sigma", "nu1", "nu2")
    spec <- astspec(data, start_pars = start_pars)
    fit <- astfit_local(spec, solver, solver_control)

    standard_errors <- sqrt(diag(solve(infoMat_ast(fit$fitted_pars))))
    fit$standard_errors <- standard_errors
    structure(fit, class = "astfit")
}

#' @rdname astfit
#' @export
summary.astfit <- function(fit) {
  fit$data <- NULL
  fit
}


astfit_local <- function(spec, solver, solver_control) {
  if (class(spec) != "astspec")
    stop("spec must be an astspec object")

  ipars <- spec$ipars
  est_idx <- which(is.na(ipars$fixed_pars))
  x0 <- ipars$start_pars[est_idx]
  names(x0) <- rownames(ipars)[est_idx]
  lb <- ipars$lower_bound[est_idx]
  ub <- ipars$upper_bound[est_idx]

  arglist <- list(data = spec$data,
                  ipars = ipars)

  # will be expanded by number of solvers developed should also implement a fixed parameter version for this
  if (solver == "nloptr") {
    res <- nloptr::nloptr(x0 = x0,
                          eval_f = llast,
                          # eval_grad_f = llast_grad,
                          lb = lb,
                          ub = ub,
                          opts = solver_control,
                          arglist = arglist)
  }

  # this is temporary, fitted should be a list with its own elements
  sol_res <- res  #list(res$pars, ...)
  fitted_pars <- res$solution
  names(fitted_pars) <- names(x0)

  return(list(data = spec$data, sol_res = sol_res, solver = solver, solver_control = solver_control,
              start_pars = spec$start_pars, fixed_pars = spec$fixed_pars, fitted_pars = fitted_pars))
}

objective_value <- function(nus, spec, solver, solver_control) {
  data <- spec$data
  ipars <- spec$ipars
  est_idx <- which(is.na(ipars$fixed_pars))
  start_pars <- ipars$start_pars[est_idx]
  names(start_pars) <- rownames(ipars)[est_idx]
  fixed_pars <- ipars$start_pars[-est_idx]
  names(fixed_pars) <- rownames(ipars)[-est_idx]

  if (nus[1] == 0 || nus[2] == 0) {
    obj = 10^5
  } else {
    spec_l <- astspec(data, start_pars = start_pars, fixed_pars = c(fixed_pars[!(names(fixed_pars) %in% c("nu1", "nu2"))], "nu1" = nus[1], "nu2" = nus[2]))
    obj <- astfit_local(spec_l, solver, solver_control)$sol_res$objective
  }
  obj
}
