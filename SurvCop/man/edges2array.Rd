% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gaussMST-cfi.r
\name{edges2array}
\alias{edges2array}
\title{Edge matrix to a vine array}
\usage{
edges2array(d, edgemat, pcv)
}
\arguments{
\item{d}{dimension or number of variables}

\item{edgemat}{edge matrix of dimension choose(d,2) x (d+1)
  with columns node1, node2, vine tree level, given[1], ..., given[d-2]
(it is not checked if this matches a proper vine array).
Each row represents an edge of the vine.}

\item{pcv}{corresponding vector of partial correlations for edges of vine}
}
\value{
Returns a list with the following named components:
\enumerate{
 \item \code{VineA} -  vine array of dimension dxd
 \item \code{pcmat} -  matrix of partial correlations by tree
}
}
\description{
Convert a vine stored as an edge matrix (columns node1, node2, level, conditioning variables) 
to a vine array
}
\details{
WARNING: It has not been checked what happens with incorrect input.
}
\examples{
edgemat6=matrix(c(
   1,2, 1, 0,0,0,0,
   1,3, 1, 0,0,0,0,
   2,4, 1, 0,0,0,0,
   1,5, 1, 0,0,0,0,
   2,6, 1, 0,0,0,0,
   2,3, 2, 1,0,0,0,
   1,4, 2, 2,0,0,0,
   3,5, 2, 1,0,0,0,
   4,6, 2, 2,0,0,0,
   3,4, 3, 2,1,0,0,
   2,5, 3, 1,3,0,0,
   1,6, 3, 2,4,0,0,
   4,5, 4, 1,3,2,0,
   3,6, 4, 2,4,1,0,
   5,6, 5, 2,4,1,3), 15,7, byrow=TRUE)
colnames(edgemat6)=c("node1","node2","level","giv1","giv2","giv3","giv4")
pcvec=c(.6,.7,.6,.7,.5, .4,.4,.3,.3, .2,.2,.1, .1,.1, .1) 
edges2array(d=6,edgemat6,pcvec)

}
