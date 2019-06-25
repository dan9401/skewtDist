#' @title AST Specification
#' the object currently have a S3 structure which may be shifted into S4 in future development
#' @description Method for creating an AST distribution object prior to fitting.
#' @param data A univariate data object, can be ... for the AST distribution to fit on.
#' @param start_pars Numeric vector of starting parameters for the optimization algorithm.
#' @param fixed_pars Numeric vector of parameters to be kept fixed during the optimization routine.
#' @name astspec
#' @examples
#' data <- rast(1000, 0.12, 0.6, 0.7, 3, 5)
#' spec <- astspec(data)

#' @rdname astspec
#' @export
astspec <- function(data, start_pars = c("mu" = 0, "sigma" = 1, "alpha" = 0.5, "nu1" = 1, "nu2" = 1),
                    fixed_pars = c()) {
    if (!is.numeric(data))
        stop("data must be numeric")
    if (length(start_pars) != 5)
        stop("start_pars must be a numeric of length 5")
    # if (!is.list(fixed_pars)) stop('fixed_pars must be a named list')

    check_bound(start_pars)
    bounds <- data.frame(name = c("mu", "sigma", "alpha", "nu1", "nu2"),
                         lower_bound = c(-Inf, 0, 0, 0, 0),
                         upper_bound = c(Inf, Inf, 1, Inf, Inf))
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
    ipars <- merge(p_df, bounds, by = "name", all = T)
    rownames(ipars) <- ipars$name
    ipars$name <- NULL

    # if (length(fixed_pars) != 0) do.call(check_bound, as.list(fixed_pars))
    structure(list(data = data, ipars = ipars), class = "astspec")
}

