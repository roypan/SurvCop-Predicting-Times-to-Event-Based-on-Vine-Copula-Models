#' van der Waerden correlation matrix when variable in the last column can be right-censored
#'
#' @description
#' van der Waerden correlation or correlation of normal scores matrix when variable in the last column can be right-censored and the other variables are mixed continuous/discrete
#' 
#' @param udata n by d matrix with values in (0,1)
#' @param is_disc logical vector of length d, TRUE in jth position for discrete variable
#' @param censor_status n by 1 logical vector, ith position TRUE if udata[i,d] is right-censored
#' @param iprint flag for printing intermediate results if TRUE
#' 
#' @return 
#' d by d positive definite van der Waerden (normal score) correlation matrix
#'
#' @example
#' d = 5
#' censoring_rate = 0.5
#' N = 1000
#' # define 5-dimensional R-vine tree structure matrix
#' Matrix = c(5, 1, 2, 3, 4,
#' 		   0, 4, 1, 2, 3,
#' 		   0, 0, 2, 1, 3,
#' 		   0, 0, 0, 3, 1,
#' 		   0, 0, 0, 0, 1)
#' Matrix = matrix(Matrix, d, d)
#' # define R-vine pair-copula family matrix
#' family = c(0, 1, 1, 1, 1,
#' 		   0, 0, 1, 1, 1,
#' 		   0, 0, 0, 1, 1,
#' 		   0, 0, 0, 0, 1,
#' 		   0, 0, 0, 0, 0)
#' family = matrix(family, d, d)
#' par = matrix(0.6, nrow = d, ncol = d)
#' par2 = matrix(0, nrow = d, ncol = d)
#' RVM = RVineMatrix(Matrix = Matrix, family = family, par = par, par2 = par2)
#' u = RVineSim(N, RVM)
#' auxiliary = pnorm(rnorm(N, -sqrt(2) * qnorm(censoring_rate)))
#' censor_status = auxiliary < u[, d]
#' u[as.logical(censor_status), d] = auxiliary[as.logical(censor_status)]
#' is_disc = c(FALSE, TRUE, TRUE, FALSE, FALSE)
#' u[, is_disc] = round(u[, is_disc])
#' vWcormat(u, is_disc, censor_status))
#'
#' @export

vWcormat=function(udata,is_disc,censor_status, iprint=FALSE){ 
	d=ncol(udata)
	rmat=matrix(0,d,d)
	rmat[1:(d-1), 1:(d-1)] = cor_mat(udata[,-d], is_disc=is_disc[-d])
	for (i in 1:(d-1)) {
		if (is_disc[i]) {
			rho_censored = vWcor_discrete(udata[,i], udata[,d], censor_status)
		} else {
			rho_censored = vWcor_continuous(udata[,i], udata[,d], censor_status)
		}
		rmat[i,d] = rho_censored
		rmat[d,i] = rho_censored
	}
	rmat[d,d] = 1
	if (!isposdef(rmat)) {
		if(iprint) { cat("not positive definite; doing a projection\n"); print(rmat); }
		rmat = as.matrix(Matrix::nearPD(rmat, corr = TRUE)$mat)
	}
	return(rmat)
}
