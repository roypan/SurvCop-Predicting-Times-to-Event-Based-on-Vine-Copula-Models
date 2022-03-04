#' predict u-score of a survival response with mixed discrete/continuous explanatory variables at a quantile
#' 
#' @description
#' predict u-score of a survival response with mixed discrete/continuous explanatory variables at quantile p such that the conditional CDF of the survival response given mixed discrete/continuous explanatory variables is equal to p
#' 
#' @param p vector of conditional CDFs to be predicted of length m
#' @param paramvec vector of length (d-1) if all the copula families associated with the survival response variable are one-parameter families; vector of length (d-1)+ otherwise
#' @param param1_complete d-1 by d-1 matrix for the first parameter of the copula models for the predictors
#' @param param2_complete d-1 by d-1 matrix for the second parameter of the copula models for the predictors
#' @param fam d by d matrix for the copula families coded by the VineCopula package
#' @param A d by d matrix for the copula structure
#' @param data_cdf_plus n by (d-1) matrix for CDF from the right
#' @param data_cdf_minus n by (d-1) matrix for CDF from the left
#' @param is_disc vector of length d representing if the variable is discrete or not
#
#' @return the conditional u-scores (n by m matrix) of the response variable given the predictors at quantiles p
#'
#' @import VineCopula
#' @export

survival_predict = function(p, paramvec, param1_complete, param2_complete, fam, A, data_cdf_plus, data_cdf_minus, is_disc, iprint=F) {
	d = ncol(A)
	n = nrow(data_cdf_plus)
	m = length(p)	# number of conditional quantiles to predict
	
	pred_result = matrix(NA, nrow = n, ncol = m)
	
	# permute all arguments based on the order of A
	A_order = diag(A)
	A = varrayperm(A, order(A_order))	 # after permutation, A has 1:d on diagonal
	
	data_cdf_plus = data_cdf_plus[, A_order[-d]]
	data_cdf_minus = data_cdf_minus[, A_order[-d]]
	
	if (n == 1) {
		data_cdf_plus = matrix(data_cdf_plus, nrow = 1)
		data_cdf_minus = matrix(data_cdf_minus, nrow = 1)
	}
	
	is_disc = is_disc[A_order]
	
	# paramvec contains parameters linked to response by tree
	paramvec = parse_parameter(paramvec, fam[1:(d-1), d], d) # make necessary transformations of the copula parameters
	param1_vec = paramvec[1:(d-1)]
	param2_vec = paramvec[d:length(paramvec)]
	
	# convert from (d-1)x(d-1) matrices to dxd matrices to include the copula parameters associated with the response variable
	parammat1 = matrix(0, nrow = d, ncol = d)
	parammat1[1:(d-1), 1:(d-1)] = param1_complete
	parammat1[1:(d-1), d] = param1_vec
	
	parammat2 = matrix(0, nrow = d, ncol = d)
	parammat2[1:(d-1), 1:(d-1)] = param2_complete
	parammat2[1:(d-1), d] = param2_vec
	
	# M and I are the matrices on Algorithm 5 of Joe (2014)
	M = matrix(0, nrow = d, ncol = d)
	for (i in 1:d) {
		for (j in i:d) {
			M[i, j] = max(A[1:i, j])
		}
	}
	
	I = matrix(0, nrow = d, ncol = d)
	for (i in 2:d) {
		for (j in i:d) {
			if (A[i, j] < M[i, j]) {
				I[i-1, M[i, j]] = 1
			}
		}
	}
	#Inputs to individual_predict are 1xd vectors because we need to make predictions for each observation by searching the solution value; solution to one observation does not necessarily help with the prediction for other observations 
	# individual_predict : dim(plus)=dim(minus)=d-1
	# individual_predict_vec : ncol(plus)=ncol(minus)=d-1
	for (k in 1:m) {
		# to avoid confusion, make arguments plus, minus instead of data_cdf_plus, data_cdf_minus
		# vec version uvec nx1, plus nx(d-1), minus nx(d-1)
		# uvec has the quantile in the interval (0,1)
		individual_predict_fortran = function(uvec, plus, minus) {
			return(.Fortran("individualpredictvec", as.integer(n), as.integer(d), as.double(parammat1), as.double(parammat2), as.integer(fam), as.integer(A), as.integer(M), as.integer(I), as.double(uvec), as.double(plus), as.double(minus), as.logical(is_disc), as.double(p[k]), output=as.double(rep(0, n)))$output)
		}
		# Conditional quantile is in the U(0,1) scale so is always between 0 and 1
		temp = bisect_extra(individual_predict_fortran, rep(.00001,n), rep(.99999,n), plus=data_cdf_plus, minus=data_cdf_minus, debug=F)
		if(iprint) { cat("n=",n,"\n"); print(temp); }
		pred_result[, k] = temp
	}
	# end of loop over quantile levels
	
	return(pred_result)	#nxm matrix, m=length(p)
}