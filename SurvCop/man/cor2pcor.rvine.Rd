% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varray-CopulaModel.r
\name{cor2pcor.rvine}
\alias{cor2pcor.rvine}
\title{correlation matrix to partial correlation vine}
\usage{
cor2pcor.rvine(rr, A)
}
\arguments{
\item{rr}{dxd positive definite correlation matrix}

\item{A}{dxd vine array with 1:d on diagonal (only upper triangle is used)}
}
\value{
List with (i) $pctree, where partial correlations for tree \eqn{\ell} are in row \eqn{\ell}, and
 (ii) $pcmat, \eqn{\rho_{jk;S}} in position (j,k) where S is conditioning vector of edge [j,k].
}
\description{
Convert a correlation matrix to a partial correlation representation based on a given vin array
}
\examples{
rr=toeplitz(c(1,.5,.25,.125,.0625,.05))
D6=Dvinearray(6)
cor2pcor.rvine(rr,D6)  # many 0s
A6=matrix(c(1,0,0,0,0,0, 1,2,0,0,0,0, 1,2,3,0,0,0, 2,1,3,4,0,0, 
         1,3,2,4,5,0, 2,4,1,3,5,6), 6,6)  
cor2pcor.rvine(rr,A6)

}
