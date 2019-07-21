#' @title AST fit
#'
#' @description Method for fitting a AST distribution
#'
#' @param data A univariate data object, can be ... for the AST distribution to fit on.
#' @param start_pars A named numeric vector of starting parameters for the optimization algorithm.
#' @param fixed_pars A named numeric vector of parameters to be kept fixed during the optimization routine.
#' @param solver Optimizer used for fitting, one of 'nloptr', 'Rsolnp' ...
#' @param solver_control Control arguments list passed to the optimizer.
#'
#' @name astfit
#'
#' @examples
#' data <- rast(1000, 0.12, 0.6, 0.3, 3, 5)
#' solver_control <- list('algorithm' = 'NLOPT_LN_BOBYQA', 'maxeval' = 1.0e5, 'xtol_rel' = 1.0e-8)
#' # solver_control <- list(eval.max = 10^3, iter.max = 10^3)
#' fit <- astfit(data, solver = 'nloptr', solver_control = solver_control)
#' # fit <- astfit(data, solver = 'nlminb', solver_control = solver_control)

#' @rdname astfit
#' @export
# fit function for ast distribution
astfit <- function(data, start_pars = c(), fixed_pars = c(), solver = c("nlminb", "nloptr", "Rsolnp"), solver_control) {
  if (!is.numeric(data) || length(data) == 0)
    stop("Data must be a numeric vector of non-zero length.")

  # function that checks parameter bounds, needs to be edited
  #check_bound(start_pars)
  #check_bound(fixed_pars)
  solver = match.arg(solver)

  fit <- astfit_local(data, start_pars, fixed_pars, solver, solver_control)
  standard_errors <- sqrt(diag(solve(infoMat_ast(fit$fitted_pars))))
  fit$standard_errors <- standard_errors

  structure(fit, class = "astfit")
}

astfit_local <- function(data, start_pars = c(), fixed_pars = c(), solver, solver_control) {
  start_time <- Sys.time()
  par_names <- c("mu", "sigma", "alpha", "nu1", "nu2")
  # insert the initial guess problem
  # possible code for using mode as the initial guess for mu
  # hist <- hist(data, breaks = 1001), only work best when a large sample of data
  # mode <- hist$mids[which(hist$counts == max(hist$counts))]
  start_pars_default <- c(mu = 0, sigma = 1, alpha = 0.5, nu1 = 2, nu2 = 2)
  start_pars <- c(start_pars, start_pars_default[!(par_names %in% names(start_pars))])

  b_df <- data.frame(name = c("mu", "sigma", "alpha", "nu1", "nu2"),
                     lower_bound = c(-Inf, 0, 0, 0, 0),
                     upper_bound = c(Inf, Inf, 1, Inf, Inf),
                     order = 1:5)
  sp_df <- data.frame(start_pars = start_pars,
                      name = names(start_pars))
  fp_df <- data.frame(fixed_pars = fixed_pars,
                      name = names(fixed_pars))
  if (length(fp_df != 0)) {
    p_df <- merge(sp_df, fp_df, by = "name", all = T)
  } else {
    sp_df$fixed_pars = NA
    p_df = sp_df
  }
  ipars <- merge(p_df, b_df, by = "name", all = T)
  ipars <- ipars[order(ipars$order), ]

  fixed_pars <- ipars$fixed_pars
  est_idx <- which(is.na(fixed_pars))
  start_pars <- ipars$start_pars[est_idx]
  lb <- ipars$lower_bound[est_idx]
  ub <- ipars$upper_bound[est_idx]
  # arglist is an argument for llast and llast_grad
  arglist <- list(data = data,
                  fixed_pars = fixed_pars)

  # will be expanded by number of solvers developed should also implement a fixed parameter version for this
  if (solver == "nloptr") {
    res <- nloptr::nloptr(x0 = start_pars,
                          eval_f = llast,
                          # eval_grad_f = llast_grad,
                          lb = lb,
                          ub = ub,
                          opts = solver_control,
                          arglist = arglist)

    # this is temporary, fitted should be a list with its own elements
    sol_res <- res  #list(res$pars, ...)
    fitted_pars <- res$solution
    names(fitted_pars) <- ipars$name[est_idx]
    objective <- sol_res$objective
  } else if (solver == "Rsolnp") {
    res <- Rsolnp::solnp(pars = start_pars,
                         fun = llast,
                         LB= lb,
                         UB = ub,
                         control = solver_control,
                         arglist = arglist)
    sol_res <- res  #list(res$pars, ...)
    fitted_pars <- res$pars
    names(fitted_pars) <- ipars$name[est_idx]
    objective <- sol_res$values[length(sol_res$values)]
  } else if (solver == "nlminb") {
    res <- nlminb(start = start_pars,
                  objective = llast,
                  gradient = llast_grad,
                  arglist = arglist,
                  control = solver_control,
                  lower = lb,
                  upper = ub)
    sol_res <- res  #list(res$pars, ...)
    fitted_pars <- res$par
    names(fitted_pars) <- ipars$name[est_idx]
    objective <- sol_res$objective
  } else {
    NA
  }
  time_elapsed <- Sys.time() - start_time

  list(data = data, sol_res = sol_res, solver = solver, solver_control = solver_control,
       start_pars = start_pars, fixed_pars = fixed_pars, fitted_pars = fitted_pars,
       objective = objective, time_elapsed = time_elapsed)
}

# log-likelihood function of the AST distributions pars: parameter values y: data which you fit the distribution on
llast <- function(pars, arglist) {
  y <- arglist$data
  fixed_pars <- arglist$fixed_pars

  est_idx <- which(is.na(fixed_pars))
  fix_idx <- which(!(is.na(fixed_pars)))
  all_pars <- c()
  all_pars[est_idx] <- pars
  all_pars[fix_idx] <- fixed_pars[fix_idx]

  mu <- all_pars[1]
  sigma <- all_pars[2]
  alpha <- all_pars[3]
  nu1 <- all_pars[4]
  nu2 <- all_pars[5]
  T_ <- length(y)
  y1 <- y[y <= mu]
  y2 <- y[y > mu]

  logl <- -T_ * log(sigma) - 0.5 * (nu1 + 1) * sum(log(1 + ((y1 - mu)/(2 * alpha * sigma * K(nu1)))^2/nu1)) -
    0.5 * (nu2 + 1) * sum(log(1 + ((y2 - mu)/(2 * (1 - alpha) * sigma * K(nu2)))^2/nu2))
  -logl
}

llast_grad <- function(pars, arglist) {
  y <- arglist$data
  fixed_pars <- arglist$fixed_pars

  est_idx <- which(is.na(fixed_pars))
  fix_idx <- which(!(is.na(fixed_pars)))
  all_pars <- c()
  all_pars[est_idx] <- pars
  all_pars[fix_idx] <- fixed_pars[fix_idx]

  mu <- all_pars[1]
  sigma <- all_pars[2]
  alpha <- all_pars[3]
  nu1 <- all_pars[4]
  nu2 <- all_pars[5]
  T_ <- length(y)
  y1 <- y[y <= mu]
  y2 <- y[y > mu]

  g_mu <- sum((nu1 + 1) / L(all_pars, y1) / nu1 * (y1 - mu) / (2 * alpha * sigma * K(nu1))^2) +
    sum((nu2 + 1) / R(all_pars, y2) / nu2 * (y2 - mu) / (2 * (1 - alpha) * sigma * K(nu2))^2)
  g_sigma<- -T_ / sigma + (nu1 + 1) / sigma * sum(1 - 1 / L(all_pars, y1)) + (nu2 + 1) / sigma * sum(1 - 1 / R(all_pars, y2))
  g_alpha <- (nu1 + 1) / alpha * sum(1 - 1 / L(all_pars, y1)) - (nu2 + 1) / (1 - alpha) * sum(1 - 1 / R(all_pars, y2))
  g_nu1 <- sum(- log(L(all_pars, y1)) / 2 + (nu1 + 1) / 2 * D(nu1) * (L(all_pars, y1) - 1) / L(all_pars, y1))
  g_nu2 <- sum(- log(R(all_pars, y2)) / 2 + (nu2 + 1) / 2 * D(nu2) * (R(all_pars, y2) - 1) / R(all_pars, y2))
  gradient <- -c(mu = g_mu, sigma = g_sigma, alpha = g_alpha, nu1 = g_nu1, nu2 = g_nu2)

  return(gradient[est_idx])
}

