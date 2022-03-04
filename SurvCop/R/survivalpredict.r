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

survivalpredict = function(p, paramvec, param1_complete, param2_complete, fam, A, data_cdf_plus, data_cdf_minus, is_disc, iprint=F) {
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
		individual_predict_vec = function(uvec, plus, minus) {
			v_plus = matrix(0,n,d)
			v_minus = matrix(0,n,d)
			vp_plus = matrix(0,n,d)
			vp_minus = matrix(0,n,d)
			s_plus = matrix(0,n,d)
			s_minus = matrix(0,n,d)
			w_plus = matrix(0,n,d)
			w_minus = matrix(0,n,d)
			
			#data_cdf_plus = c(data_cdf_plus, u)
			#data_cdf_minus = c(data_cdf_minus, u)
			plus = cbind(plus, uvec)
			minus = cbind(minus, uvec)
			
			# initialization steps
			for (j in 1:d) {
				s_plus[,j] = plus[,A[1, j]]
				s_minus[,j] = minus[,A[1, j]]
				w_plus[,j] = plus[,j]
				w_minus[,j] = minus[,j]
			}
			
			# levels 2:(d-1)
			for (l in 2:(d-1)) {
				#s_plus[s_plus >= .999] = .999
				#s_minus[s_minus >= .999] = .998
				#w_plus[w_plus >= .999] = .999
				#w_minus[w_minus >= .999] = .998
				#s_plus[s_plus <= .001] = .002
				#s_minus[s_minus <= .001] = .001
				#w_plus[w_plus <= .001] = .002
				#w_minus[w_minus <= .001] = .001
				for (j in l:d) {
					if (I[l-1, j] == 1) {
						if (is_disc[j]) {
							vp_plus[,j] = (BiCopCDF(s_plus[,j], w_plus[,j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j]) - BiCopCDF(s_plus[,j], w_minus[,j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])) / (w_plus[,j] - w_minus[,j])
							vp_minus[,j] = (BiCopCDF(s_minus[,j], w_plus[,j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j]) - BiCopCDF(s_minus[,j], w_minus[,j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])) / (w_plus[,j] - w_minus[,j])
						} else {
							vp_plus[,j] = BiCopHfunc2(s_plus[,j], w_plus[,j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])
							vp_minus[,j] = BiCopHfunc2(s_minus[,j], w_plus[,j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])
						}
					}
					if (is_disc[A[l-1, j]]) {
						v_plus[,j] = (BiCopCDF(s_plus[,j], w_plus[,j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j]) - BiCopCDF(s_minus[,j], w_plus[,j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])) / (s_plus[,j] - s_minus[,j])
						v_minus[,j] = (BiCopCDF(s_plus[,j], w_minus[,j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j]) - BiCopCDF(s_minus[,j], w_minus[,j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])) / (s_plus[,j] - s_minus[,j])
					} else {
						v_plus[,j] = BiCopHfunc2(w_plus[,j], s_plus[,j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])
						v_minus[,j] = BiCopHfunc2(w_minus[,j], s_plus[,j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])
					}
				}
				
				for (j in (l+1):d) {
					if (A[l, j] == M[l, j]) {
						s_plus[,j] = v_plus[,M[l, j]]
						s_minus[,j] = v_minus[,M[l, j]]
					} else {
						s_plus[,j] = vp_plus[,M[l, j]]
						s_minus[,j] = vp_minus[,M[l, j]]
					}
					w_plus[,j] = v_plus[,j]
					w_minus[,j] = v_minus[,j]
				}
				
				if (l == (d - 1)) {
					if (is_disc[A[d-1, d]]) {
						ucond = (BiCopCDF(s_plus[,d], w_plus[,d], fam[d-1, d], parammat1[d-1, d], parammat2[d-1, d]) - BiCopCDF(s_minus[,d], w_plus[,d], fam[d-1, d], parammat1[d-1, d], parammat2[d-1, d])) / (s_plus[,d] - s_minus[,d])
					} else {
						ucond = BiCopHfunc2(w_plus[,d], s_plus[,d], fam[d-1, d], parammat1[d-1, d], parammat2[d-1, d])
					}
				}
			}
			# ucond and result are nx1 vectors
			result = ucond - p[k]
			#if (is.na(result)) {
			#	if (u <= .1) result = -1
			#	if (u >= .9) result = 1
			#}
			return(result)
		}
		# end of function individual_predict_vec 
		
		# Conditional quantile is in the U(0,1) scale so is always between 0 and 1
		temp = bisect_extra(individual_predict_vec, rep(.00001,n), rep(.99999,n), plus=data_cdf_plus, minus=data_cdf_minus, debug=F)
		if(iprint) { cat("n=",n,"\n"); print(temp); }
		pred_result[, k] = temp
	}
	# end of loop over quantile levels
	
	return(pred_result)	#nxm matrix, m=length(p)
}