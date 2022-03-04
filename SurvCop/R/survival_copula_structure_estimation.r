#' vine copula structure estimation with a survival response
#' 
#' @description
#' this function estimates the vine structure (such that the response variable is contained in a leaf node at every level of the vine), selects the bivariate copula families, and estimates the parameters associated with the selected copulas for a survival response variable with mixed discrete/continuous explanatory variables
#' 
#' @param data_cdf: n by d matrix for u-scores
#' @param censor_status: vector of length n representing if the response variable is censored or not
#' @param is_disc: vector of length d representing if the variable is discrete or not
#' @param prlevel: print level for nlm
#' @param all_fam: all bivariate copula families to be considered as coded in VineCopula package
#
#' @return a list with the following attributes:
#' \enumerate{
#'  \item \code{best_param} - estimated copula parameters associated with the censored response variable given the best vine structure and family matrix found
#'  \item \code{param1_complete} - matrix of the first parameters associated with the copula modesl for the predictors 
#'  \item \code{param2_complete} - matrix of the second parameters associated with the copula modesl for the predictors 
#'  \item \code{fam_estimated} - family matrix for the copula models
#'  \item \code{A_mst} - vine structure matrix
#' }
#'
#' @import VineCopula
#' @export

survival_copula_structure_estimation = function(data_cdf, censor_status, is_disc, prlevel = 2, all_fam = c(1, 2, 4, 14, 5, 7, 10, 20, 17)) {
	d = ncol(data_cdf)
	n = nrow(data_cdf)
	
	data_cdf_plus = data_cdf
	data_cdf_minus = data_cdf
	for (j in 1:d) {
		if (is_disc[j]) {
			data_cdf_unique = c(0, sort(unique(data_cdf[, j])))
			for (i in 1:n) {
				data_cdf_minus[i, j] = data_cdf_unique[which(data_cdf_unique == data_cdf[i, j]) - 1]
			}
		}
	}
	
	# Find the vine structure using Bo's method with the correlation matrix estimated assuming bivariate Gaussian copula everywhere
	cormat = matrix(0, d, d)
	cormat[1:(d-1), 1:(d-1)] = cor_mat(data_cdf[, -d], is_disc = is_disc[-d])

	for (i in 1:(d-1)) {
		if (is_disc[i]) {
			rho_censored = vWcor_discrete(data_cdf[, i], data_cdf[, d], censor_status)
		} else {
			rho_censored = vWcor_continuous(data_cdf[, i], data_cdf[, d], censor_status)
		}
		cormat[i, d] = rho_censored
		cormat[d, i] = rho_censored
	}
	cormat[d, d] = 1

	# check positive definiteness
	if (!isposdef(cormat)) {
		cormat = as.matrix(Matrix::nearPD(cormat, corr = TRUE)$mat) # check if it is true with some correlation matrices
	}
	print(cormat)
	
	A_mst = get_vine_array(cormat, n = n, cfi_bd = 1.01)$vine_array
	
	# use the VineCopula package to select the copula models on each edge for the given vine structure using all the fully observed variables
	RVM_complete = rvine_cop_select_mix(data_cdf[, 1:(d-1)], is_disc[1:(d-1)], familyset = c(1,2,4,5,7,14,17,10,20), A_mst[(d-1):1, (d-1):1])
	print(RVM_complete)

	param1_complete = RVM_complete$par[(d-1):1, (d-1):1]
	param2_complete = RVM_complete$par2[(d-1):1, (d-1):1]
	fam_complete = RVM_complete$family[(d-1):1, (d-1):1]
	Matrix_complete = RVM_complete$Matrix[(d-1):1, (d-1):1]
	
	# test all Gaussian baseline vine copula model
	fam = matrix(0, nrow = d, ncol = d)
	fam[1:(d-1), 1:(d-1)] = fam_complete
	fam[1:(d-1), d] = 1

	modelfit = nlm(f = survival_nllk, p = rep(0, d-1), param1_complete = param1_complete, param2_complete = param2_complete, fam = fam, A = A_mst, data_cdf_plus = data_cdf_plus, data_cdf_minus = data_cdf_minus, censor_status = censor_status, is_disc = is_disc, stepmax = 10, iterlim = 300, print.level=prlevel)

	param_est = modelfit$estimate
	rho_est = parse_parameter(param_est, rep(1, d - 1), d)[1:(d-1)]

	# if |parcial correlation| > .3, test alternative copula families
	partial_cor = cor2pcor.rvine(cormat, A_mst)$pctree[1:(d-1), d]
	partial_cor_ind = abs(partial_cor) > .3

	indices = rep(1, d-1)
	aic_min = rep(Inf, length(all_fam))
	best_param = param_est
	param_est_list = list()

	# sequentially choose the best family to use on each edge
	for (i in 1:(d - 1)) {
		if (abs(partial_cor_ind[i]) > .3 ) {
			alt_indices = matrix(indices, nrow = length(all_fam), ncol = d-1, byrow = TRUE)
			alt_indices[, i] = all_fam
			for (j in 1:length(all_fam)) {
				starting_val = best_param
				fam = matrix(0, nrow = d, ncol = d)
				fam[1:(d-1), 1:(d-1)] = fam_complete
				fam[1:(d-1), d] = alt_indices[j, ]
				for (k in i:(d - 1)) {
					if (alt_indices[j, k] == 1) { #Gaussian
						#pass
					} else if (alt_indices[j, k] == 2) { # t
						starting_val = c(best_param, 5)
					} else if (alt_indices[j, k] == 4) { # gumbel
						starting_val[k] = depmeas2cpar(rho_est[i], type = "rhoN", copname = "gumbel")
					} else if (alt_indices[j, k] == 14) { # reflected gumbel
						starting_val[k] = depmeas2cpar(rho_est[i], type = "rhoN", copname = "gumbel")
					} else if (alt_indices[j, k] == 5) { # frank
						starting_val[k] = depmeas2cpar(rho_est[i], type = "rhoN", copname = "frank")
					} else if (alt_indices[j, k] == 7) { # bb1 bb1r
						bb1_param = bb1.tau2eqlm(bvn.cpar2tau(tanh(starting_val[i] / 10)))
						starting_val[k] = bb1_param[1]
						starting_val = c(starting_val, bb1_param[2])
					}
				}
				modelfit = nlm(f = survival_nllk2, p = starting_val, param1_complete = param1_complete, param2_complete = param2_complete, fam = fam, A = A_mst, data_cdf_plus = data_cdf_plus, data_cdf_minus = data_cdf_minus, censor_status = censor_status, is_disc = is_disc, stepmax = 5, iterlim = 300, print.level=prlevel)
				aic_min[j] = 2*length(starting_val) + 2*modelfit$minimum
				param_est_list[[j]] = modelfit$estimate
			}
			best_param = param_est_list[[which.min(aic_min)]]
			indices[i] = all_fam[which.min(aic_min)]
		}
	}

	fam_estimated = fam
	fam_estimated[1:(d-1), d] = indices
	print(indices)
	
	best_param[indices %in% c(1,2)] = tanh(best_param[indices %in% c(1,2)] / 10)
	
	return_list = list('best_param' = best_param, 'param1_complete' = param1_complete, 'param2_complete' = param2_complete, 'fam_estimated' = fam_estimated, 'A_mst' = A_mst)
	return(return_list)
}
