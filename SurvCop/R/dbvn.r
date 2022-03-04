#' bivariate normal density (internal function)
#'
#' @description
#' bivariate normal density with zero means and unit variances
#'
#' @param x vector of length 2 or 2-column matrix
#' @param rho correlation parameter with -1 < rho < 1
#' @param ilog TRUE for log density, FALSE for density
#'
#' @return bivariate standard normal density with correlation rho, or log density 
#'
#' @examples
#' rho = .8
#' N = 1000
#' x = qnorm(VineCopula::BiCopSim(N, 1, rho, par2 = 0, obj = NULL, check.pars = TRUE))
#' dbvn(x, rho)

dbvn=function(x, rho, ilog=FALSE) { 
	if(is.matrix(x)) { x1=x[,1]; x2=x[,2] }
	else { x1=x[1]; x2=x[2] }
	qf=x1^2+x2^2-2*rho*x1*x2
	qf=qf/(1-rho^2)
	con=sqrt(1-rho^2)*(2*pi)
	if(!ilog) pdf=exp(-0.5*qf)/con
	else pdf=-0.5*qf-log(con)
	return(pdf)
}
