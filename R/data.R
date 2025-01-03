#' Weekly return data for Large and Small Caps stocks from 1997-01-01 to 2010-12-31.
#'
#' @format
#' \describe{
#'   \item{retLW}{A xts object of 730 rows and 20 columns, each column represents a different large-cap stock.}
#'   \item{retSW}{A xts object of 730 rows and 20 columns, each column represents a different small-cap stock.}
#' }
#'
#' @source \url{https://finance.yahoo.com/lookup/}
#'
#' @name weeklyRtnData
#' @aliases retLW, retSW
#' @usage data(retLW)
#'        data(retSW)
#' @importFrom xts xts
"_PACKAGE"

#' @rdname weeklyRtnData
"retLW"

#' @rdname weeklyRtnData
"retSW"