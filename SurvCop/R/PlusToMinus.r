#' Plus to Minus for discrete variables
#' 
#' @description
#' Given a discrete Uscore (or cdf) matrix F^+_i = P(X_i <= x_i), 
#' calculate the left limit of the cdf F^-_i = P(X_i < x_i).
#' Internal function.
#'
#' @param dataDisc A discrete Uscore matrix.
#' @param listSupport A list of sorted support points, that is the possible values 
#' each variable can take. If not input, the function will produce it.
#' 
#' @return Returns a list with the following named components:
#' \enumerate{
#' \item \code{dataMinus} - A discrete Uscore (cdf) matrix wih the same dimension as dataDisc
#' \item \code{listSupport} - List of support points for each discrete variable
#' }
#'
PlusToMinus <- function (dataDisc, listSupport = NULL) {
  if (length(dataDisc) == 0)
    return(dataDisc)
  
  if (is.vector(dataDisc))
    dataDisc <- matrix(dataDisc, ncol = 1)
  
  if (!is.matrix(dataDisc))
    dataDisc <- as.matrix(dataDisc)
  
  if (is.null(listSupport)) {
    listDataDisc <- split(dataDisc, col(dataDisc))
    listSupport <- lapply(listDataDisc, FUN = function(x) sort(unique(x)))
  }
  
  d <- ncol(dataDisc)
  stopifnot(is.list(listSupport))
  stopifnot(d == length(listSupport))
  
  result <- dataDisc
  for (k in 1:d) {
    tryCatch(
      result[, k] <-
        unlist(lapply(dataDisc[, k], 
                      FUN = function (x) {
                        ind <- which.min(abs(x - listSupport[[k]]))
                        ifelse(ind > 1, 
                               listSupport[[k]][ind - 1],
                               0)
                      })),
      error = function (e) stop("listSupport does not match input data.")
    )
  }
  
  return(list(
    dataMinus = result,
    listSupport = listSupport))
}
