#' truncated vine based on comparative fit index
#' 
#' @description
#' Apply CFI (comparative fit index) calculations to a correlation matrix and a vine array to get a truncated vine
#' 
#' @param rmat positive definite dxd correlation matrix, with d>2 
#' @param A  dxd vine array (only upper triangle is used)
#' @param n  sample size used to compute rmat
#' @param CFIbd lower bound of CFI (comparative fit index) to truncate, default 0.95, set as 1.01 for no truncation
#' @param iprint  set as T for intermediate prints, default=F
#' 
#' @return Returns a list with the following named components:
#' \enumerate{
#'  \item \code{VineA} -  vine array of dimension dxd
#'  \item \code{pcmat} -  matrix of partial correlations by tree
#'  \item \code{treeweight} -  vector of length d-1 with
#'                 \eqn{\sum_e \log(1-\rho_e^2)} for trees 1,...d-1
#'  \item \code{edgemat} -  matrix with columns node1 node2 level vector-of-conditioning
#'  \item \code{fitval} - vector with cumsum(treeweight)/sum(treeweight)
#'  \item \code{CFIv} - vector of CFI values
#'  \item \code{ntrunc} - truncation level to reach CFIbd
#' }
#' 
#' @examples
#' rmat=matrix(c(
#'   1.000,  0.172, -0.062,  0.385,  0.499, 0.267,  0.578,
#'   0.172,  1.000, -0.047,  0.493,  0.602, 0.396,  0.600,
#'  -0.062, -0.062,  1.000, -0.120, -0.092, 0.070, -0.070,
#'   0.385,  0.493, -0.120,  1.000,  0.607, 0.557,  0.742,
#'   0.499,  0.602, -0.092,  0.607,  1.000, 0.672,  0.883,
#'   0.267,  0.396,  0.070,  0.557,  0.672, 1.000,  0.685,
#'   0.578,  0.600, -0.070,  0.742,  0.883, 0.685,  1.000), 7,7)
#' A7=matrix(c( 
#'   5, 5, 5, 5, 5, 4, 5,
#'   0, 4, 4, 4, 4, 5, 4,
#'   0, 0, 1, 1, 6, 6, 1,
#'   0, 0, 0, 6, 1, 1, 6,
#'   0, 0, 0, 0, 2, 2, 2,
#'   0, 0, 0, 0, 0, 3, 3,
#'   0, 0, 0, 0, 0, 0, 7), 7,7, byrow=TRUE)
#' trvine.CFI(rmat, A7, n=400, CFIbd=0.98)
#' 
#' @export
#' 
trvine.CFI=function(rmat, A, n, CFIbd, iprint=FALSE)
{ d=nrow(rmat)
  if(d<=2) return(rmat)
  perm=diag(A); rr=rmat[perm,perm]
  iperm=order(perm)
  AA=varrayperm(A,iperm)
  #print(AA)  # should be 1:d on diagonal
  #print(rr)
  pcobj=cor2pcor.rvine(rr,AA); pc=pcobj$pctree
  #print(pc)
  logdet = log(det(rr))
  #logdet2 = sum(log(1-pc^2))
  #cat(logdet2,logdet,"\n")  # should be same to a few decimal places
  treeweight = rep(0,d-1)
  CFIv=rep(0,d-1)
  nu0=(d*(d-1))/2
  D0=-n*logdet
  sumtreewt=0
  for(itree in 1:(d-1))
  { treeweight[itree] = sum(log(1-pc[itree,(itree+1):d]^2))
    sumtreewt = sumtreewt+treeweight[itree]
    Dt=n*(sumtreewt-logdet); nut=((d-itree)*(d-itree-1))/2
    numer=max(0,Dt-nut); denom=max(numer,D0-nu0)
    CFIv[itree]=1-numer/denom
    if(iprint) cat("CFI", itree, Dt,nut, numer, denom, 1-numer/denom,"\n")
    if(CFIv[itree]>=CFIbd) break
  }
  ntrunc=itree
  fitval=cumsum(treeweight)
  fitval=fitval/logdet
  Atrunc=A
  if(ntrunc<d-1) 
  { for(i in (itree+1):(d-1)) Atrunc[i,(i+1):d]=0 }
  list(VineA=Atrunc,pcmat=pc[1:ntrunc,],
     treeweight=treeweight[1:ntrunc],fitval=fitval[1:ntrunc], 
     CFIv=CFIv[1:ntrunc],ntrunc=ntrunc)
}

#' Get CFI for a partial correlation array based on a truncation level
#'
#' @description
#' A helper function that computes the comparative fit index (CFI)
#' for a given correlation matrix and partial correlation array and truncation level
#'
#' @param rmat Correlation matrix.
#' @param pc Partial correlation array (computed using get_cfi_from_vine_array)
#' @param n Sample size.
#' @param trunc_level Truncation level.
#'
#' @return CFI
#'
#' @details
#' CFI (comparative fit index) is defined as follows.
#' CFI = 1- numerator/denominator,\cr
#' numerator = \eqn{\max(0,D_t-\nu_t)},\cr
#' denominator = \eqn{\max(0,D_t-\nu_t,D_0-\nu_0)},\cr
#' \eqn{D_0=-n\log(\det(R))},\cr
#' \eqn{D_t = n[-L_t(V) - \log(\det(R))]} is decreasing as t increases,\cr
#' \eqn{L_t(V) = -\sum_{1:t} \sum_e -\log(1-r_e^2)} (sum pcor up to tree t),\cr
#' \eqn{\nu_t = (d-t)(d-t-1)/2; \nu_0=d(d-1)/2},\cr
#' \eqn{L_t(V)} is increasing in t.
#' 
get_cfi_from_pc <- function(rmat, pc, n, trunc_level) {
  d <- nrow(rmat)
  D_0 <- -n * log(det(rmat))
  nu_0 <- d * (d - 1) / 2
  pc <- pc[1:trunc_level, ]
  L_t <- sum(-log(1 - pc[!is.na(pc) & pc != 0] ^ 2))
  D_t <- n * (-log(det(rmat)) - L_t)
  nu_t <- (d - trunc_level) * (d - trunc_level - 1) / 2
  1 - max(0, D_t - nu_t) / max(0, D_0 - nu_0, D_t - nu_t)
}

#' Get CFI from a vine array based on a truncation level
#'
#' @description
#' A helper function that computes the comparative fit index (CFI)
#' for a given correlation matrix and vine array and truncation level
#'
#' @param rmat Correlation matrix.
#' @param vine_array Vine array.
#' @param n Sample size.
#' @param trunc_level Truncation level.
#'
#' @return CFI
#'
get_cfi_from_vine_array <- function(rmat, vine_array, n, trunc_level) {
  diag_varray <- diag(vine_array)
  perm <- order(diag_varray)
  vine_array_permed <- varrayperm(vine_array, perm)
  pctree <- cor2pcor.rvine(rmat[diag_varray, diag_varray], vine_array_permed)$pctree
  get_cfi_from_pc(rmat, pctree, n, trunc_level)
}
