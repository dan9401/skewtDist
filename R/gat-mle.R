#' @title Fitting function for Asymmetric Student-t distribution
#'
#' @name gatMLE
#' @aliases gatmle
#'
#' @description Method for fitting an gat distribution to a univariate data series by Maximum Likelihood Estimation,
#' returns an \code{gat} object.
#'
#' @param data a univariate data object to be fitted
#' @param start_pars a named numeric vector of starting parameters for the optimization algorithm, not all parameters are needed
#' @param fixed_pars a named numeric vector of parameters to be kept fixed during the optimization routine, not all parameters are needed
#' @param solver solver used for MLE, one of 'nlminb', 'nloptr', 'Rsolnp', default is 'nlminb'
#' @param solver_control list of control arguments passed to the solver
#' @param symmetric a logical argument, when TRUE, the function fits an SST distribution(Symmetric Student-t, nu1 = nu2) instead, default to FALSE
#'
#' @return
#' A \code{gat} object(S3), the components of the object are:
#'     \item{data}{the univariate data object for the gat distribution to be fitted}
#'     \item{solver}{the solver called}
#'     \item{solver_control}{the list of control argumetns passed to the solver called}
#'     \item{start_pars}{named numeric vector of starting parameters used}
#'     \item{fixed_pars}{named numeric vector of fixed parameters used}
#'     \item{symmetric}{logical argument controlling the symmetry of tail parameters in the MLE}
#'     \item{solver_result}{output of the called solver}
#'     \item{fitted_pars}{named vector of fitted arguemnts of the gat distribution}
#'     \item{objective}{the optimal log-likelihood value obtained by the solver}
#'     \item{time_elapsed}{the time elapesed for the MLE routine}
#'     \item{message}{the message of convergence status produced by the called solver}
#'     \item{standard_errors}{standard errors of the fitted parameters}
#'
#' @details
#' The \code{gatMLE} function fits an gat distribution to a univariate data series by estimating the distribution parameters
#' through Maximum Likelihood Estimation.
#'
#' For details of the list of control arguments, please refer to \code{nlminb}, \code{nloptr::nloptr}, \code{Rsolnp::solnp}
#'
#' @references
#' Zhu, D., & Galbraith, J. W. (2010). A generalized asymmetric Student-t distribution with application to financial econometrics. Journal of Econometrics, 157(2), 297-305.\url{https://www.sciencedirect.com/science/article/pii/S0304407610000266}
#' \url{https://econpapers.repec.org/paper/circirwor/2009s-13.htm}
#'
#' @examples
#' pars <- c(0, 1, 1.5, 1.2, 2, 4)
#' data <- rgat(1000, pars = pars)
#' fit <- gatMLE(data)
#' fit
#' 
#' @importFrom stats nlminb

#' @rdname gatMLE
#' @export
gatMLE <- function(data, start_pars = c(), fixed_pars = c(), solver = c("nlminb", "nloptr", "Rsolnp"), solver_control = list()) {
  if (!is.numeric(data) || length(data) == 0)
    stop("Data must be a numeric vector of non-zero length.")

  # function that checks parameter bounds, needs to be edited
  #check_bound(start_pars)
  #check_bound(fixed_pars)
  solver = match.arg(solver)

  fit <- gatfit_local(data, start_pars, fixed_pars, solver, solver_control)
  standard_errors <- sqrt(diag(solve(gatInfoMat(fit$fitted_pars, data = data)))/length(data))
  fit$standard_errors <- standard_errors

  structure(fit, class = "gat")
}

gatfit_local <- function(data, start_pars, fixed_pars, solver, solver_control) {
  start_time <- Sys.time()
  par_names <- c("mu", "phi", "alpha", "r", "c", "nu")
  # insert the initial guess problem
  # possible code for using mode as the initial guess for mu
  # hist <- hist(data, breaks = 1001), only work best when a large sample of data
  # mode <- hist$mids[which(hist$counts == max(hidocumentst$counts))]
  start_pars_default <- c(mu = 0, phi = 1, alpha = 1, r = 1, c = 1, nu = 2)
  start_pars <- c(start_pars, start_pars_default[!(par_names %in% names(start_pars))])

  b_df <- data.frame(name = c("mu", "phi", "alpha", "r", "c", "nu"),
                     lower_bound = c(-Inf, 0, 0, 0, 0, 0),
                     upper_bound = c(Inf, Inf, Inf, Inf, Inf, Inf),
                     order = 1:6)
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
  # arglist is an argument for llgat and llgat_grad
  arglist <- list(data = data,
                  fixed_pars = fixed_pars)

  # will be expanded by number of solvers developed should also implement a fixed parameter version for this
  if (solver == "nloptr") {
    res <- nloptr::nloptr(x0 = x0,
                          eval_f = llgat,
                          # eval_grad_f = llgat_grad,
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
                         fun = llgat,
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
                  objective = llgat,
                  gradient = llgat_grad,
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

  list(data = data, solver = solver, solver_control = solver_control,
       start_pars = start_pars, fixed_pars = fixed_pars,
       solver_result = sol_res, fitted_pars = fitted_pars,
       objective = objective, time_elapsed = time_elapsed, message )
}



llgat <- function(pars, arglist) {
  y <- arglist$data
  fixed_pars <- arglist$fixed_pars

  est_idx <- which(is.na(fixed_pars))
  fix_idx <- which(!(is.na(fixed_pars)))
  all_pars <- c()
  all_pars[est_idx] <- pars
  all_pars[fix_idx] <- fixed_pars[fix_idx]

  mu <- all_pars[1]
  phi <- all_pars[2]
  alpha <- all_pars[3]
  r <- all_pars[4]
  c <- all_pars[5]
  nu <- all_pars[6]
  T_ <- length(y)

  z <- (y - mu) / phi
  g <- z + sqrt(1 + z^2)
  p <- nu / (alpha * (1 + r^2))
  q <- p * r^2
  A <- (c*g)^(alpha*r) + (c*g)^(-alpha/r)

  logl <- T_ * log( alpha * (1 + r^2) / (r * phi) ) - sum( nu / alpha * log(A) ) - T_ * log(beta(p, q)) - sum( 1/2 * log(1 + z^2) )
  -logl
}

llgat_grad <- function(pars, arglist) {
  y <- arglist$data
  fixed_pars <- arglist$fixed_pars

  est_idx <- which(is.na(fixed_pars))
  fix_idx <- which(!(is.na(fixed_pars)))
  all_pars <- c()
  all_pars[est_idx] <- pars
  all_pars[fix_idx] <- fixed_pars[fix_idx]

  mu <- all_pars[1]
  phi <- all_pars[2]
  alpha <- all_pars[3]
  r <- all_pars[4]
  c <- all_pars[5]
  nu <- all_pars[6]
  T_ <- length(y)

  z <- (y - mu) / phi
  g <- z + sqrt(1 + z^2)
  p <- nu / (alpha * (1 + r^2))
  q <- p * r^2
  A <- (c*g)^(alpha*r) + (c*g)^(-alpha/r)
  dAdg <- alpha*r*(c*g)^(alpha*r)/g - alpha/r*(c*g)^(-alpha/r)/g
  dAdr <- (c*g)^(alpha*r)*alpha*log(c*g) + (c*g)^(-alpha/r)*alpha*log(c*g)/r^2
  dgdz <- 1 + z/sqrt(1+z^2)
  dzdmu <- -1/phi
  dzdphi <- -(y-mu)/phi^2
  dpdalpha <- -p/alpha
  dpdnu <- p/nu
  dpdr <- -2*r/(1+r^2)*p
  dqdalpha <- -q/alpha
  dqdnu <- q/nu
  dqdr <- 2/(r*(1+r^2))*q
  dbdp <- beta(p,q)*(digamma(p) - digamma(p+q))
  dbdq <- beta(p,q)*(digamma(q) - digamma(p+q))
  dbdnu <- dbdq * dqdnu + dbdp * dpdnu
  dbdalpha <- dbdq * dqdalpha + dbdp * dpdalpha
  dbdr <- dbdq * dqdr + dbdp * dpdr

  g_mu <- sum( -nu/alpha/A*dAdg*dgdz*dzdmu - z/(1+z^2)*dzdmu )
  g_phi<- -T_/phi + sum( -nu/alpha/A*dAdg*dgdz*dzdphi - z/(1+z^2)*dzdphi )
  g_alpha <-  T_/alpha + sum(nu/alpha^2*log(A) - nu/alpha/A*( r*(c*g)^(alpha*r)*log(c*g) - 1/r*(c*g)^(-alpha/r)*log(c*g) ) ) - T_/beta(p, q) * dbdalpha
  g_r <- T_*2*r/(1+r^2) - T_/r + sum( -nu/alpha/A*dAdr ) - T_/beta(p,q) * dbdr
  g_c <- sum( -nu/alpha/A* (alpha*r*g^(alpha*r)*c^(alpha*r-1) - alpha/r*g^(-alpha/r)*c^(-alpha/r-1)) )
  g_nu <- sum( -1/alpha*log(A) ) - T_/beta(p,q) * dbdnu
  gradient <- -c(mu = g_mu, phi = g_phi, alpha = g_alpha, g_r = g_r, g_c = g_c, nu = g_nu)

  return(gradient[est_idx])
}
