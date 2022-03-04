#' conditional CDF of the survival response given mixed discrete/continuous explanatory variables
#' 
#' @description
#' conditional CDF of the survival response given mixed discrete/continuous explanatory variables
#' 
#' @param u-score scalar or vector of the response variable of length n
#' @param paramvec vector of length (d-1) if all the copula families associated with the survival response variable are one-parameter families; vector of length (d-1)+ otherwise
#' @param param1_complete d-1 by d-1 matrix for the first parameter of the copula models for the predictors
#' @param param2_complete d-1 by d-1 matrix for the second parameter of the copula models for the predictors
#' @param fam d by d matrix for the copula families coded by the VineCopula package
#' @param A d by d matrix for the copula structure
#' @param data_cdf_plus n by (d-1) matrix for CDF from the right
#' @param data_cdf_minus n by (d-1) matrix for CDF from the left
#' @param is_disc vector of length d representing if the variable is discrete or not
#
#' @return the conditional CDF of the response variable given the explanatory variables of length n
#'
#' @import VineCopula
#' @export

survivalconditionalcdf = function(u, paramvec, param1_complete, param2_complete, fam, A, data_cdf_plus, data_cdf_minus, is_disc) {
	d = ncol(A)
	n = nrow(data_cdf_plus)
	
	# permute all arguments based on the order of A
	A_order = diag(A)
	A = varrayperm(A, order(A_order))
	
	data_cdf_plus = data_cdf_plus[, A_order[-d]]
	data_cdf_minus = data_cdf_minus[, A_order[-d]]
	
	if (n == 1) {
		data_cdf_plus = matrix(data_cdf_plus, nrow = 1)
		data_cdf_minus = matrix(data_cdf_minus, nrow = 1)
	}
	
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
	
	v_plus = rep(0, d)
	v_minus = rep(0, d)
	vp_plus = rep(0, d)
	vp_minus = rep(0, d)
	s_plus = rep(0, d)
	s_minus = rep(0, d)
	w_plus = rep(0, d)
	w_minus = rep(0, d)
	
	result = numeric(length(u))
	
	for (i in 1:length(u)) {
		data_cdf_plus_i = c(data_cdf_plus[i, ], u[i])
		data_cdf_minus_i = c(data_cdf_minus[i, ], u[i])
		
		# initialization steps
		for (j in 1:d) {
			s_plus[j] = data_cdf_plus_i[A[1, j]]
			s_minus[j] = data_cdf_minus_i[A[1, j]]
			w_plus[j] = data_cdf_plus_i[j]
			w_minus[j] = data_cdf_minus_i[j]
		}
		
		if (d == 2) {
			if (is_disc[A[d-1, d]]) {
				ucond = (BiCopCDF(s_plus[d], w_plus[d], fam[d-1, d], parammat1[d-1, d], parammat2[d-1, d]) - BiCopCDF(s_minus[d], w_plus[d], fam[d-1, d], parammat1[d-1, d], parammat2[d-1, d])) / (s_plus[d] - s_minus[d])
			} else {
				# condition on first variable (not discrete)
				ucond = BiCopHfunc2(w_plus[d], s_plus[d], fam[d-1, d], parammat1[d-1, d], parammat2[d-1, d])
			}
		} else {
			# levels 2:(d-1)
			for (l in 2:(d-1)) {
				s_plus[s_plus >= 1] = .999
				s_minus[s_minus >= .999] = .998
				w_plus[w_plus >= 1] = .999
				w_minus[w_minus >= .999] = .998
				s_plus[s_plus <= .001] = .002
				s_minus[s_minus <= 0] = .001
				w_plus[w_plus <= .001] = .002
				w_minus[w_minus <= 0] = .001
				for (j in l:d) {
					if (I[l-1, j] == 1) {
						if (is_disc[j]) {
							vp_plus[j] = (BiCopCDF(s_plus[j], w_plus[j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j]) - BiCopCDF(s_plus[j], w_minus[j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])) / (w_plus[j] - w_minus[j])
							vp_minus[j] = (BiCopCDF(s_minus[j], w_plus[j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j]) - BiCopCDF(s_minus[j], w_minus[j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])) / (w_plus[j] - w_minus[j])
						} else {
							# condition on second variable (not discrete) for vprime 
							vp_plus[j] = BiCopHfunc2(s_plus[j], w_plus[j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])
							vp_minus[j] = BiCopHfunc2(s_minus[j], w_plus[j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])
						}
					}
					if (is_disc[A[l-1, j]]) {
						v_plus[j] = (BiCopCDF(s_plus[j], w_plus[j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j]) - BiCopCDF(s_minus[j], w_plus[j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])) / (s_plus[j] - s_minus[j])
						v_minus[j] = (BiCopCDF(s_plus[j], w_minus[j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j]) - BiCopCDF(s_minus[j], w_minus[j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])) / (s_plus[j] - s_minus[j])
					} else {
						# condition on first variable (not discrete)
						v_plus[j] = BiCopHfunc2(w_plus[j], s_plus[j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])
						v_minus[j] = BiCopHfunc2(w_minus[j], s_plus[j], fam[l-1, j], parammat1[l-1, j], parammat2[l-1, j])
					}
				}
				
				for (j in (l+1):d) {
					if (A[l, j] == M[l, j]) {
						s_plus[j] = v_plus[M[l, j]]
						s_minus[j] = v_minus[M[l, j]]
					} else {
						s_plus[j] = vp_plus[M[l, j]]
						s_minus[j] = vp_minus[M[l, j]]
					}
					w_plus[j] = v_plus[j]
					w_minus[j] = v_minus[j]
				}
				
				if (l == (d - 1)) {
					if (is_disc[A[d-1, d]]) {
						ucond = (BiCopCDF(s_plus[d], w_plus[d], fam[d-1, d], parammat1[d-1, d], parammat2[d-1, d]) - BiCopCDF(s_minus[d], w_plus[d], fam[d-1, d], parammat1[d-1, d], parammat2[d-1, d])) / (s_plus[d] - s_minus[d])
					} else {
						# condition on first variable (not discrete)
						ucond = BiCopHfunc2(w_plus[d], s_plus[d], fam[d-1, d], parammat1[d-1, d], parammat2[d-1, d])
					}
				}
			}
		}
		result[i] = ucond
	}
	return(result)
}
