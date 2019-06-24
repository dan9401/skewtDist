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

    ipars <- spec$ipars
    est_idx <- which(is.na(ipars$fixed_pars))
    x0 <- ipars$start_pars[est_idx]
    names(x0) <- rownames(ipars)[est_idx]
    lb <- ipars$lower_bound[est_idx]
    ub <- ipars$upper_bound[est_idx]

    # insert grid search here, maybe

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

    structure(list(data = spec$data, start_pars = spec$start_pars, fixed_pars = spec$fixed_pars, fitted_pars = fitted_pars,
        solver = solver, solver_control = solver_control, sol_res = sol_res), class = "astfit")
}
