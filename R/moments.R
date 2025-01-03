#' Calculate all moments of a distribution
#' 
#' @name moments
#' 
#' @description The generic S3 methods for calculating all moments of a fitted distribution.
#' 
#' @param x A AST/GAT fit object of class \code{\link{ast}} / \code{\link{gat}}.
#' @param ... Additional arguments passed into the calculating functions of moments.

#' @rdname moments
#' @export
moments <- function(x, ...) {
  UseMethod("moments", x)
}