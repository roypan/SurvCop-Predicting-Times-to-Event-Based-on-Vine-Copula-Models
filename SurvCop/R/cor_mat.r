#' correlation matrix : vander Waarden/ polychoric/ polyserial
#' 
#' @description
#' A correlation matrix of mixed continuous/discrete data.
#' For continuous-continuous, calculate the vander Waarden correlation of the normal scores.
#' For continuous-discrete, calculate the polyserial correlation using \code{polycor::polychor}.
#' For discrete-discrete, calculate the polychoric correlation using \code{polychor::polyserial}.
#'
#' @param data dataframe or matrix
#' @param is_disc logical vector, TRUE for the i-th variable if it is discrete. By default all FALSE.
#' @param debug debug flag
#'
#' @return correlation matrix
#'
#' @examples
#' x = matrix(rnorm(700), nrow = 100, ncol = 7)
#' x[, c(1,2,3)] = round(x[, c(1,2,3)])
#' cor_mat(x, c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
#'
#' @import polycor
#' @export
#'
cor_mat <- function(data, is_disc = rep(FALSE, ncol(data)), debug = FALSE) {
  d <- ncol(data)
  stopifnot(d == length(is_disc))
  rmat <- matrix(0, d, d)
  
  nscore_continuous <- nscore(data[, !is_disc])
  rmat[!is_disc, !is_disc] <- cor(nscore_continuous)
  
  for (i in which(is_disc)) {
    for (j in 1:d) {
      if (i == j)
        next
      
      if (!is_disc[j]) {
        # i is discrete and j is continuous, polyserial
        rmat[i, j] <-
          rmat[j, i] <- polycor::polyserial(data[, i], data[, j])
      } else if (is_disc[j]) {
        # i and j are both discrete, polychoric
        rmat[i, j] <-
          rmat[j, i] <- polycor::polychor(data[, i], data[, j])
      }
    }
  }
  
  # If rmat is not positive definite, use Matrix::nearPD to find the nearest
  # positive definite matrix
  if (!isposdef(rmat))
    rmat <- as.matrix(Matrix::nearPD(rmat)$mat)
  diag(rmat) <- 1
  
  if (debug) {
    cat("\n\nCorrelation matrix:\n")
    print(rmat)
  }
  
  return(rmat)
}

# From CopulaModel
#' Is positive definite?
#' 
#' @description
#' Check if a square symmetric matrix is positive definite.
#'
#' @param amat symmetric matrix
#'
#' @return TRUE if amat is positive definite, FALSE otherwise
#'
#' @examples
#' a1=matrix(c(1,.5,.5,1),2,2); print(isposdef(a1))
#' a2=matrix(c(1,1.5,1.5,1),2,2); print(isposdef(a2))
#'
#' @export
#'
isposdef = function(amat)
{
  tt = try(chol(amat), silent = T)
  ifelse("matrix" %in% class(tt), T, F)
}
