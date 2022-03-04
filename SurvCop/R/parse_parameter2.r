#' parse bivariate copula parameters for likelihood estimation (internal function)
#' 
#' @description
#' Parse bivariate copula parameters for likelihood estimation so that all parameters are within the range. Note that parameters for Gaussian and t copulas before parsing are in the range of (-infinity, infinity), parameters for Gaussian and t copulas after parsing are in the range of (-1, 1). A tanh(x/10) transformation is taken.
#' 
#' @param paramvec vector of copula parameters to be parsed, the first d-1 elements are the first parameters for the copula models, the remaining elements (if any) are the second parameters for the copula models
#' @param fam_vec vector of copula families of length d-1
#' @param d dimension of the dataset
#' 
#' @return the parsed parameter vector

parse_parameter2 = function(paramvec, fam_vec, d) {
	fam_vec_twoparam = fam_vec %in% c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)
	param1_vec = paramvec[1:(d-1)]
	param2_vec = rep(0, d-1)
	if (length(paramvec) > (d-1)) {
		param2_vec[fam_vec_twoparam] = paramvec[d:length(paramvec)]
	}
	for (i in 1:(d-1)) {
		if (fam_vec[i] == 1) { # Gaussian
			# tanh transformation
			param1_vec[i] = tanh(param1_vec[i]/10)
			if (param1_vec[i] <= -1) param1_vec[i] = -.999
			if (param1_vec[i] >= 1) param1_vec[i] = .999
		} else if (fam_vec[i] == 2) { # t
			# tanh transformation
			param1_vec[i] = tanh(param1_vec[i]/10)
			if (param1_vec[i] <= -1) param1_vec[i] = -.999
			if (param1_vec[i] >= 1) param1_vec[i] = .999
			if (param2_vec[i] <= 2) param2_vec[i] = 2.001
		} else if (fam_vec[i] %in% c(3, 13, 23, 33)) { # Clayton
			if (param1_vec[i] <= 0) param1_vec[i] = 0.001
			if (param1_vec[i] > 100) param1_vec[i] = 100
		} else if (fam_vec[i] %in% c(4, 14, 24, 34)) { # Gumbel
			if (param1_vec[i] < 1) param1_vec[i] = 1.001
			if (param1_vec[i] > 20) param1_vec[i] = 20
		} else if (fam_vec[i] == 5) { # Frank
			if (param1_vec[i] < -30) param1_vec[i] = -30
			if (param1_vec[i] > 30) param1_vec[i] = 30
		} else if (fam_vec[i] %in% c(6, 16, 26, 36)) { # Joe
			if (param1_vec[i] <= 1) param1_vec[i] = 1.001
			if (param1_vec[i] > 50) param1_vec[i] = 50
		} else if (fam_vec[i] %in% c(7, 17, 27, 37)) { # BB1
			if (param1_vec[i] <= 0) param1_vec[i] = .001
			if (param1_vec[i] > 7) param1_vec[i] = 7
			if (param2_vec[i] < 1) param2_vec[i] = 1
			if (param2_vec[i] > 7) param2_vec[i] = 7
		} else if (fam_vec[i] %in% c(8, 18, 28, 38)) { # BB6
			if (param1_vec[i] < 1) param1_vec[i] = 1
			if (param1_vec[i] > 6) param1_vec[i] = 6
			if (param2_vec[i] < 1) param2_vec[i] = 1
			if (param2_vec[i] > 8) param2_vec[i] = 8
		} else if (fam_vec[i] %in% c(9, 19, 29, 39)) { # BB7
			if (param1_vec[i] < 1) param1_vec[i] = 1
			if (param1_vec[i] > 6) param1_vec[i] = 6
			if (param2_vec[i] <= 0) param2_vec[i] = 0.001
			if (param2_vec[i] > 75) param2_vec[i] = 75
		} else if (fam_vec[i] %in% c(10, 20, 30, 40)) { # BB8
			if (param1_vec[i] < 1) param1_vec[i] = 1
			if (param1_vec[i] > 8) param1_vec[i] = 8
			if (param2_vec[i] < 1e-4) param2_vec[i] = 1e-4
			if (param2_vec[i] > 1) param2_vec[i] = 1
		}
	}
	return(c(param1_vec, param2_vec))
}