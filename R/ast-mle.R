#' @title Fitting function for Asymmetric Student-t distribution
#'
#' @name astMLE
#' @aliases astmle
#'
#' @description Method for fitting an AST distribution to a univariate data series by Maximum Likelihood Estimation,
#' returns an \code{ast} object.
#'
#' @param data a univariate data object to be fitted
#' @param start_pars a named numeric vector of starting parameters for the optimization algorithm, not all parameters are needed
#' @param fixed_pars a named numeric vector of parameters to be kept fixed during the optimization routine, not all parameters are needed
#' @param solver solver used for MLE, one of 'nlminb', 'nloptr', 'Rsolnp', default is 'nlminb'
#' @param solver_control list of control arguments passed to the solver
#' @param symmetric a logical argument, when TRUE, the function fits an SST distribution(Symmetric Student-t, nu1 = nu2) instead, default to FALSE
#'
#' @return
#' A \code{ast} object(S3), the components of the object are:
#'     \item{data}{the univariate data object for the AST distribution to be fitted}
#'     \item{solver}{the solver called}
#'     \item{solver_control}{the list of control argumetns passed to the solver called}
#'     \item{start_pars}{named numeric vector of starting parameters used}
#'     \item{fixed_pars}{named numeric vector of fixed parameters used}
#'     \item{symmetric}{logical argument controlling the symmetry of tail parameters in the MLE}
#'     \item{solver_result}{output of the called solver}
#'     \item{fitted_pars}{named vector of fitted arguemnts of the AST distribution}
#'     \item{objective}{the optimal log-likelihood value obtained by the solver}
#'     \item{time_elapsed}{the time elapesed for the MLE routine}
#'     \item{message}{the message of convergence status produced by the called solver}
#'     \item{standard_errors}{standard errors of the fitted parameters}
#'
#' @details
#' The \code{astMLE} function fits an AST distribution to a univariate data series by estimating the distribution parameters
#' through Maximum Likelihood Estimation.
#'
#' For details of the list of control arguments, please refer to \code{nlminb}, \code{nloptr::nloptr}, \code{Rsolnp::solnp}
#'
#' @references
#' Zhu, D., & Galbraith, J. W. (2010). A generalized asymmetric Student-t distribution with application to financial econometrics. Journal of Econometrics, 157(2), 297-305.\url{https://www.sciencedirect.com/science/article/pii/S0304407610000266}
#' \url{https://econpapers.repec.org/paper/circirwor/2009s-13.htm}
#'
#' @examples
#' pars <- c(0.12, 0.6, 0.6, 3, 5)
#' data <- rast(1000, pars = pars)
#' 
#' solver_control <- list(eval.max = 10^3, iter.max = 10^3)
#' fit <- astMLE(data, solver_control = solver_control)
#' fit
#' 
#' @importFrom stats nlminb

#' @rdname astMLE
#' @export
# fit function for ast distribution
astMLE <- function(data, start_pars = c(), fixed_pars = c(), solver = c("nlminb", "nloptr", "Rsolnp"), solver_control = list(), symmetric = FALSE) {
  if (!is.numeric(data) || length(data) == 0)
    stop("Data must be a numeric vector of non-zero length.")

  # function that checks parameter bounds, needs to be edited
  #check_bound(start_pars)
  #check_bound(fixed_pars)
  solver = match.arg(solver)
  fit <- astFit_local(data, start_pars, fixed_pars, solver, solver_control, symmetric)
  if (symmetric == TRUE) {
    standard_errors <- sqrt(diag(solve(sstInfoMat(fit$fitted_pars)))/length(data))
  } else {
    standard_errors <- sqrt(diag(solve(astInfoMat(fit$fitted_pars)))/length(data))
  }

  fit$standard_errors <- standard_errors

  structure(fit, class = "ast")
}

astFit_local <- function(data, start_pars, fixed_pars, solver, solver_control, symmetric) {
  start_time <- Sys.time()
  ipars <- i_pars(start_pars, fixed_pars, symmetric)

  fixed_pars <- ipars$fixed_pars
  if (is.null(start_pars)) {
    start_pars <- rep(NA, 5)
  } else {
    start_pars <- start_pars[ipars$name]
  }
  names(fixed_pars) <- names(start_pars) <- ipars$name
  est_idx <- which(is.na(fixed_pars))
  x0 <- ipars$start_pars[est_idx]
  lb <- ipars$lower_bound[est_idx]
  ub <- ipars$upper_bound[est_idx]
  # arglist is an argument for llast and llast_grad
  arglist <- list(data = data,
                  fixed_pars = fixed_pars)

  # will be expanded by number of solvers developed should also implement a fixed parameter version for this
  if (solver == "nloptr") {
    res <- nloptr::nloptr(x0 = x0,
                          eval_f = llast,
                          eval_grad_f = llast_grad,
                          lb = lb,
                          ub = ub,
                          opts = solver_control,
                          arglist = arglist)
    # this is temporary, fitted should be a list with its own elements
    sol_res <- res  #list(res$pars, ...)
    fitted_pars <- res$solution
    names(fitted_pars) <- ipars$name[est_idx]
    objective <- sol_res$objective
    message <- sol_res$message
  } else if (solver == "Rsolnp") {
    res <- Rsolnp::solnp(pars = x0,
                         fun = llast,
                         LB= lb,
                         UB = ub,
                         control = solver_control,
                         arglist = arglist)
    sol_res <- res  #list(res$pars, ...)
    fitted_pars <- res$pars
    names(fitted_pars) <- ipars$name[est_idx]
    objective <- sol_res$values[length(sol_res$values)]
    message <- sol_res$convergence
  } else if (solver == "nlminb") {
    res <- nlminb(start = x0,
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
    message <- sol_res$message
  } else {
    NA
  }
  time_elapsed <- Sys.time() - start_time

  fitted_pars <- c(fitted_pars, fixed_pars[!is.na(fixed_pars)])
  if (symmetric == TRUE) {
    fitted_pars <- fitted_pars[c("mu", "sigma", "alpha", "nu")]
  } else {
    fitted_pars <- fitted_pars[c("mu", "sigma", "alpha", "nu1", "nu2")]
  }

  list(data = data, solver = solver, solver_control = solver_control,
       start_pars = start_pars, fixed_pars = fixed_pars, symmetric = symmetric,
       solver_result = sol_res, fitted_pars = fitted_pars,
       objective = objective, time_elapsed = time_elapsed, message = message)
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
  nu2 <- ifelse(length(all_pars) == 4, nu1, all_pars[5])
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
  nu2 <- ifelse(length(all_pars) == 4, nu1, all_pars[5])
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

  # gradient <- pmax(gradient, -1e6)
  # gradient <- pmin(gradient, 1e6)
  return(gradient[est_idx])
}

i_pars <- function(start_pars, fixed_pars, symmetric) {

  # insert the initial guess problem
  # possible code for using mode as the initial guess for mu
  # hist <- hist(data, breaks = 1001), only work best when a large sample of data
  # mode <- hist$mids[which(hist$counts == max(hist$counts))]
  if (symmetric == TRUE) {
    b_df <- data.frame(name = c("mu", "sigma", "alpha", "nu"),
                       lower_bound = c(-1000, 0, 0, 0),
                       upper_bound = c(1000, 1000, 1, 1000),
                       order = 1:4)
    start_pars_default <- c(mu = 0, sigma = 1, alpha = 0.5, nu = 2)
  } else {
    b_df <- data.frame(name = c("mu", "sigma", "alpha", "nu1", "nu2"),
                       lower_bound = c(-Inf, 0, 0, 0, 0),
                       upper_bound = c(Inf, Inf, 1, Inf, Inf),
                       order = 1:5)
    start_pars_default <- c(mu = 0, sigma = 1, alpha = 0.5, nu1 = 2, nu2 = 2)
  }
  start_pars <- c(start_pars, start_pars_default[!(b_df$name %in% names(start_pars))])
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
  # return
  ipars
}
