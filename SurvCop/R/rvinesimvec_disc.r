#' simulation for R-vine copulas with mixed discrete/continuous variables
#' 
#' @description
#' simulating u-scores from a R-vine copulas with mixed discrete/continuous variables based on rosenblatt transformation
#' 
#' @param nsim sample size for simulation
#' @param A d*d vine array, or ntrunc*d vine array as only ntrunc rows are used
#' @param fam d by d matrix for the copula families coded by the VineCopula package
#' @param param1 d by d matrix for the first parameter of the copula models
#' @param param2 d by d matrix for the second parameter of the copula models
#' @param is_disc vector of length d representing if the variable is discrete or not
#' @param iprint print flag for intermediate results
#' @param discrete_values possible values for the discrete variables to take at the original scale
#' @param marginal_dist string of the pdf function of the marginal distribution (for discrete variables)
#' @param marginal_param vector of the parameters of the marginal distribution (for discrete variables)
#' 
#' @return the simulated u-scores
#' 
#' @export

rvinesimvec_disc = function(nsim, A, fam, param1, param2, is_disc, discrete_values = c(-1.5, -.5, .5, 1.5), marginal_dist = 'pnorm', marginal_param = c(0, 1), iprint = F) {
	d = ncol(A)
	pdist = get(marginal_dist)
	
	# transform the discrete values to cut points
	cut_points = c(-Inf)
	for (i in 1:(length(discrete_values) - 1)) {
		cut_points = c(cut_points, (discrete_values[i] + discrete_values[i+1]) / 2)
	}
	if (length(marginal_param) == 1){
		cut_points = pdist(c(cut_points, Inf), marginal_param)
	} else if (length(marginal_param) == 2){
		cut_points = pdist(c(cut_points, Inf), marginal_param[1], marginal_param[2])
	} else if (length(marginal_param) == 3){
		cut_points = pdist(c(cut_points, Inf), marginal_param[1], marginal_param[2], marginal_param[3])
	}
	print(cut_points)
	
	u = matrix(runif(nsim * d), nrow = nsim, ncol = d)
	
	which_j = A[1, 1]
	if (!is_disc[which_j]) {
		data_cdf_plus = matrix(u[, which_j], ncol = 1)
		data_cdf_minus = matrix(u[, which_j], ncol = 1)
	} else {
		u[, which_j] = cut(u[, which_j], cut_points, labels = FALSE)
		emp_cdf = ecdf(u[, which_j])
		data_cdf_plus = matrix(pmin(emp_cdf(u[, which_j]), .999), ncol = 1)
		data_cdf_minus = pmax(matrix(emp_cdf(u[, which_j] - 1), ncol = 1), 0)
		u[, which_j] = data_cdf_plus
	}
	
	for (j in 2:d) {
		A_sub = A[1:j, 1:j]
		which_j = A[j, j]
		diagA_sub = diag(A_sub)
		dict = data.frame(Col1=c(0, diagA_sub), Col2=0:j)
		A_sub = matrix(dict$Col2[match(A_sub, dict$Col1)], nrow = j, byrow = FALSE)
		param1_j = param1[1:(j-1), j]
		#if (sum(fam[1:(j-1), j] %in% c(1, 2)) > 0) param1_j[fam[1:(j-1), j] %in% c(1, 2)] = atanh(param1_j[fam[1:(j-1), j] %in% c(1, 2)]) * 10
		if (!is_disc[which_j]) {
			for (i in 1:nsim) {
				u[i, which_j] = survival_predict(u[i, which_j], paramvec = c(param1_j, param2[1:(j-1), j][param2[1:(j-1), j] != 0]), param1_complete = param1[1:(j-1), 1:(j-1)], param2_complete = param2[1:(j-1), 1:(j-1)], fam = fam[1:j, 1:j], A = A_sub, data_cdf_plus = matrix(data_cdf_plus[i, ], nrow = 1), data_cdf_minus = matrix(data_cdf_minus[i, ], nrow = 1), is_disc = is_disc[diagA_sub])
			}
			data_cdf_plus = cbind(data_cdf_plus, u[, which_j])
			data_cdf_minus = cbind(data_cdf_minus, u[, which_j])
		} else {
			conditional_cdf = matrix(0, nrow = nsim, ncol = length(cut_points))
			conditional_cdf[, length(cut_points)] = 1
			data_cdf_plus_j = c()
			data_cdf_minus_j = c()
			for (i in 1:nsim) {
				for (k in 2:(length(cut_points) - 1)) {
					conditional_cdf[i, k] = survival_conditional_cdf(cut_points[k], paramvec = c(param1_j, param2[1:(j-1), j][param2[1:(j-1), j] != 0]), param1_complete = param1[1:(j-1), 1:(j-1)], param2_complete = param2[1:(j-1), 1:(j-1)], fam = fam[1:j, 1:j], A = A_sub, data_cdf_plus = matrix(data_cdf_plus[i, ], nrow = 1), data_cdf_minus = matrix(data_cdf_minus[i, ], nrow = 1), is_disc = c(is_disc[diagA_sub[1:(j-1)]], FALSE))
				}
				#print(conditional_cdf[i, ])
				ind = max(which(u[i, which_j] >= conditional_cdf[i, ]))
				data_cdf_plus_j = c(data_cdf_plus_j, ind)
				data_cdf_minus_j = c(data_cdf_minus_j, ind - 1)
			}
			emp_cdf = ecdf(data_cdf_plus_j)
			data_cdf_plus = cbind(data_cdf_plus, pmin(emp_cdf(data_cdf_plus_j), .999))
			data_cdf_minus = cbind(data_cdf_minus, pmax(emp_cdf(data_cdf_minus_j), 0))
			u[, which_j] = pmin(emp_cdf(data_cdf_plus_j), .999)
		}
	}
	return(u)
}