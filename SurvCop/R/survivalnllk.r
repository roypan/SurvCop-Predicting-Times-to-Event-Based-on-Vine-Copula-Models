#' negative log-likelihood of vine copula model for survival response
#' 
#' @description
#' negative log-likelihood of a vine copula model for a survival response (the last column of the data matrix) and mixed discrete/continuous explanatory variables (the other columns of the data matrix)
#' 
#' @param paramvec vector of length (d-1) if all the copula families associated with the survival response variable are one-parameter families; vector of length (d-1)+ otherwise
#' @param param1_complete d-1 by d-1 matrix for the first parameter of the copula models for the predictors
#' @param param2_complete d-1 by d-1 matrix for the second parameter of the copula models for the predictors
#' @param fam d by d matrix for the copula families coded by the VineCopula package
#' @param A d by d matrix for the copula structure
#' @param data_cdf_plus n by (d-1) matrix for CDF from the right
#' @param data_cdf_minus n by (d-1) matrix for CDF from the left
#' @param censor_status vector of length n representing if the response variable is censored or not for each observation
#' @param is_disc vector of length d representing if each variable is discrete or not
#'
#' @return the negative log-likelihood of the specified vine copula model for a survival response
#'
#' @import VineCopula
#' @export

survivalnllk = function(paramvec, param1_complete, param2_complete, fam, A, data_cdf_plus, data_cdf_minus, censor_status, is_disc) {
	d = ncol(A)
	n = nrow(data_cdf_plus)
	
	# permute all arguments based on the order of A
	A_order = diag(A)
	A = varrayperm(A, order(A_order))
	
	data_cdf_plus = data_cdf_plus[, A_order]
	data_cdf_minus = data_cdf_minus[, A_order]
	is_disc = is_disc[A_order]
	
	paramvec = parse_parameter(paramvec, fam[1:(d-1), d], d)
	param1_vec = paramvec[1:(d-1)]
	param2_vec = paramvec[d:length(paramvec)]
	
	parammat1 = matrix(0, nrow = d, ncol = d)
	parammat1[1:(d-1), 1:(d-1)] = param1_complete
	parammat1[1:(d-1), d] = param1_vec
	
	parammat2 = matrix(0, nrow = d, ncol = d)
	parammat2[1:(d-1), 1:(d-1)] = param2_complete
	parammat2[1:(d-1), d] = param2_vec
	
	v_plus = matrix(0, nrow = n, ncol = d)
	v_minus = matrix(0, nrow = n, ncol = d)
	vp_plus = matrix(0, nrow = n, ncol = d)
	vp_minus = matrix(0, nrow = n, ncol = d)
	s_plus = matrix(0, nrow = n, ncol = d)
	s_minus = matrix(0, nrow = n, ncol = d)
	w_plus = matrix(0, nrow = n, ncol = d)
	w_minus = matrix(0, nrow = n, ncol = d)
	
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
	
	nllk = 0
	
	# initialization steps
	for (j in 1:d) {
		s_plus[, j] = data_cdf_plus[, A[1, j]]
		s_minus[, j] = data_cdf_minus[, A[1, j]]
		w_plus[, j] = data_cdf_plus[, j]
		w_minus[, j] = data_cdf_minus[, j]
	}
	
	# first level
	# uncensored
	for (j in 2:(d-1)) {
		nllk = nllk - sum(mixed_copula_llk(s_plus[, j], s_minus[, j], w_plus[, j], w_minus[, j], is_disc[A[1, j]], is_disc[j], fam[1, j], parammat1[1, j], parammat2[1, j]))
	}
	#print(nllk)
	# censored
	nllk = nllk - sum(mixed_copula_llk(s_plus[!censor_status, d], s_minus[!censor_status, d], w_plus[!censor_status, d], w_minus[!censor_status, d], is_disc[A[1, d]], is_disc[d], fam[1, d], parammat1[1, d], parammat2[1, d]))
	#print(nllk)
	
	# levels 2:(d-1)
	for (l in 2:(d-1)) {
		s_plus[s_plus >= .999] = .999
		s_minus[s_minus >= .999] = .998
		w_plus[w_plus >= .999] = .999
		w_minus[w_minus >= .999] = .998
		s_plus[s_plus <= .001] = .002
		s_minus[s_minus <= .001] = .001
		w_plus[w_plus <= .001] = .002
		w_minus[w_minus <= .001] = .001
		for (j in l:d) {
			if (I[l-1, j] == 1) {
				if (is_disc[j]) {
					vp_plus[, j] = (BiCopCDF(s_plus[, j], w_plus[, j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j]) - BiCopCDF(s_plus[, j], w_minus[, j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])) / (w_plus[, j] - w_minus[, j])
					vp_minus[, j] = (BiCopCDF(s_minus[, j], w_plus[, j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j]) - BiCopCDF(s_minus[, j], w_minus[, j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])) / (w_plus[, j] - w_minus[, j])
				} else {
					# condition on second variable for vp=vprime
					vp_plus[, j] = BiCopHfunc2(s_plus[, j], w_plus[, j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])
					vp_minus[, j] = BiCopHfunc2(s_minus[, j], w_plus[, j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])
				}
			}
			if (is_disc[A[l-1, j]]) {
				v_plus[, j] = (BiCopCDF(s_plus[, j], w_plus[, j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j]) - BiCopCDF(s_minus[, j], w_plus[, j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])) / (s_plus[, j] - s_minus[, j])
				v_minus[, j] = (BiCopCDF(s_plus[, j], w_minus[, j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j]) - BiCopCDF(s_minus[, j], w_minus[, j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])) / (s_plus[, j] - s_minus[, j])
			} else {
				# condition on first variable for v
				v_plus[, j] = BiCopHfunc2(w_plus[, j], s_plus[, j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])
				v_minus[, j] = BiCopHfunc2(w_minus[, j], s_plus[, j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])
			}
		}
		
		for (j in (l+1):d) {
			if (A[l, j] == M[l, j]) {
				s_plus[, j] = v_plus[, M[l, j]]
				s_minus[, j] = v_minus[, M[l, j]]
			} else {
				s_plus[, j] = vp_plus[, M[l, j]]
				s_minus[, j] = vp_minus[, M[l, j]]
			}
			w_plus[, j] = v_plus[, j]
			w_minus[, j] = v_minus[, j]
		}
		
		if (l < (d - 1)) {
			# uncensored
			for (j in (l+1):(d-1)) {
				nllk = nllk - sum(mixed_copula_llk(s_plus[, j], s_minus[, j], w_plus[, j], w_minus[, j], is_disc[A[l, j]], is_disc[j], fam[l, j], parammat1[l, j], parammat2[l, j]))
			}
			#print(nllk)
			# censored
			nllk = nllk - sum(mixed_copula_llk(s_plus[!censor_status, d], s_minus[!censor_status, d], w_plus[!censor_status, d], w_minus[!censor_status, d], is_disc[A[l, d]], is_disc[d], fam[l, d], parammat1[l, d], parammat2[l, d]))
			#print(nllk)
		} else {
			# uncensored
			nllk = nllk - sum(mixed_copula_llk(s_plus[!censor_status, d], s_minus[!censor_status, d], w_plus[!censor_status, d], w_minus[!censor_status, d], is_disc[A[d-1, d]], is_disc[d], fam[d-1, d], parammat1[d-1, d], parammat2[d-1, d]))
			#print(nllk)
			if (is_disc[A[d-1, d]] & sum(censor_status) > 0) {
				# censored
				nllk = nllk - sum(log(1 - (BiCopCDF(s_plus[censor_status, d], w_plus[censor_status, d], fam[d-1, d], parammat1[d-1, d], parammat2[d-1, d]) - BiCopCDF(s_minus[censor_status, d], w_plus[censor_status, d], fam[d-1, d], parammat1[d-1, d], parammat2[d-1, d])) / (s_plus[censor_status, d] - s_minus[censor_status, d])))
				#print(nllk)
			} else {
				# censored
				nllk = nllk - sum(log(1 - BiCopHfunc2(w_plus[censor_status, d], s_plus[censor_status, d], fam[d-1, d], parammat1[d-1, d], parammat2[d-1, d])))
				#print(nllk)
			}
		}
	}
	return(nllk)
}