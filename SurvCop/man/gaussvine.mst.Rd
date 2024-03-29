% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gaussMST-cfi.r
\name{gaussvine.mst}
\alias{gaussvine.mst}
\title{Sequential MST based on partial correlations}
\usage{
gaussvine.mst(rmat, n, CFIbd = 0.95, iprint = FALSE, iprint2 = FALSE)
}
\arguments{
\item{rmat}{correlation matrix}

\item{n}{sample size}

\item{CFIbd}{lower bound of CFI (comparative fit index) to stop, default 0.95, set as 1.01 for no truncation}

\item{iprint}{print flag for intermediate results on truncated vine}

\item{iprint2}{print flag for intermediate results on graph objects}
}
\value{
Returns a list with the following named components:
\enumerate{
 \item \code{VineA} -  vine array of dimension dxd
 \item \code{pcmat} -  matrix of partial correlations by tree
 \item \code{treeweight} -  vector of length d-1 with
                sum_edge log(1-rho[edge]^2) for trees 1,...d-1
 \item \code{edgemat} -  matrix with columns node1 node2 level vector-of-conditioning
 \item \code{fitval} - cumsum(treeweight)/\eqn{sum_{1:(d-1)}} treeweight
 \item \code{CFIv} - vector of CFI values
 \item \code{ntrunc} - truncation level to reach CFIbd
}
}
\description{
Sequential MST based on partial correlations for best truncated Gaussian vine
This depends on the R package igraph
for the minimum spanning tree algorithm and operations on graphs.
}
\details{
CFI = 1- numerator/denominator, \cr
numerator = \eqn{\max(0,D_t-\nu_t)},\cr
denominator = \eqn{\max(0,D_t-\nu_t,D_0-\nu_0)},\cr
\eqn{D_0=-n\log(\det(R))},\cr
\eqn{D_t = n[-L_t(V) - \log(\det(R))]} is decreasing as t increases,\cr
\eqn{L_t(V) = -\sum_{1:t} \sum_e -\log(1-r_e^2)} (sum pcor up to tree t),\cr
\eqn{\nu_t = (d-t)(d-t-1)/2; \nu_0=d(d-1)/2},\cr
\eqn{L_t(V)} is increasing in t.
}
\examples{
n=400
rmat=matrix(c(
  1.000,  0.172, -0.062,  0.385,  0.499, 0.267,  0.578,
  0.172,  1.000, -0.047,  0.493,  0.602, 0.396,  0.600,
 -0.062, -0.062,  1.000, -0.120, -0.092, 0.070, -0.070,
  0.385,  0.493, -0.120,  1.000,  0.607, 0.557,  0.742,
  0.499,  0.602, -0.092,  0.607,  1.000, 0.672,  0.883,
  0.267,  0.396,  0.070,  0.557,  0.672, 1.000,  0.685,
  0.578,  0.600, -0.070,  0.742,  0.883, 0.685,  1.000),7,7)
vinestr095 = gaussvine.mst(rmat, n, CFIbd=0.95, iprint=TRUE)
print(vinestr095$VineA)
print(vinestr095$CFIv)
vinestr099 = gaussvine.mst(rmat, n, CFIbd=0.99, iprint=TRUE)
print(vinestr099$VineA)
print(vinestr099$CFIv)
vinestr100 = gaussvine.mst(rmat, n, CFIbd=1.01, iprint=FALSE)
print(vinestr100$VineA)
print(vinestr100$CFIv)

}
\references{
For use of CFI in the context of partial correlation vines, please see: \cr
Brechmann E C and Joe H (2015), Truncation of vine copulas using fit indices.
Journal of Multivariate Analysis, 138, 19-33.
}
