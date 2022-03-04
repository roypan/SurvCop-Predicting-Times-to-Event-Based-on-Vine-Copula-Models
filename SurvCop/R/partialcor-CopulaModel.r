# functions for partial correlations, originally from CopulaModel library
# allpcor   .. not needed
# trvine.logdet .. not needed
# subset2index .. not needed
# d2brev .. not needed



#============================================================

# S = covariance or correlation matrix
# given = vector indices for the given or conditioning variables
# j,k = indices for the conditioned variables
# Output: partial correlation of variables j,k given indices in 'given'


#' Partial correlation
#'
#' @description
#' Computation of a partial correlation from a covariance matrix
#'
#' @param S covariance or correlation matrix, of dimension dxd
#' @param given vector indices for the given or conditioning variables
#' @param j,k indices for the conditioned variables
#'
#' @return partial correlation of variables j,k given indices in 'given'
#'
#' @details
#' WARNING: It is not checked that j,k,given are indices in 1:d, where
#'   d is the dimension. It is not checked that intersect(c(j,k),given) is empty
#'
#' @examples
#' rr=toeplitz(c(1,.5,.25,.125,.05))
#' partcor(rr,c(1),3,4)
#' partcor(rr,c(1),3,5)
#' partcor(rr,c(1,2),3,4)
#' partcor(rr,c(1,2),3,5)
#'
#' @export
#'
partcor = function(S,given,j,k)
{
  S11 = S[given,given]
  jk = c(j,k)
  S12 = S[given,jk]
  S21 = S[jk,given]
  S22 = S[jk,jk]
  if (length(given) > 1) {
    tem = solve(S11,S12); Om212 = S21 %*% tem
  }
  else {
    tem = S12 / S11; Om212 = outer(S21,tem)
  }
  om11 = 1 - Om212[1,1]
  om22 = 1 - Om212[2,2]
  om12 = S[j,k] - Om212[1,2]
  return( om12 / sqrt(om11 * om22) )
}


