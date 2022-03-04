#' log-likelihood for bivariate copula with mixed discrete/continuous variables
#' 
#' @description
#' log-likelihood for bivariate copula with mixed discrete/continuous variables
#' 
#' @param u_plus scalar or vector, F+ (u-score/CDF from right side) for the first variable
#' @param u_minus scalar or vector, F- (u-score/CDF from left side) for the first variable
#' @param v_plus scalar or vector, F+ (u-score/CDF from right side) for the second variable
#' @param v_minus scalar or vector, F- (u-score/CDF from left side) for the second variable
#' @param is_disc_u logical value indicating if the first variable is discrete
#' @param is_disc_v logical value indicating if the second variable is discrete
#' @param fam scalar representing bivariate copula family that is coded 
#' families: \cr
#' \code{0} = independence copula \cr
#' \code{1} = Gaussian copula \cr
#' \code{2} = Student t copula (t-copula) \cr
#' \code{3} = Clayton copula \cr
#' \code{4} = Gumbel copula \cr
#' \code{5} = Frank copula \cr
#' \code{6} = Joe copula \cr
#' \code{7} = BB1 copula \cr
#' \code{8} = BB6 copula \cr
#' \code{9} = BB7 copula \cr
#' \code{10} = BB8 copula \cr
#' \code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' \code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' \code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' \code{17} = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
#' \code{18} = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
#' \code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
#' \code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
#' \code{23} = rotated Clayton copula (90 degrees) \cr
#' \code{24} = rotated Gumbel copula (90 degrees) \cr
#' \code{26} = rotated Joe copula (90 degrees) \cr
#' \code{27} = rotated BB1 copula (90 degrees) \cr
#' \code{28} = rotated BB6 copula (90 degrees) \cr
#' \code{29} = rotated BB7 copula (90 degrees) \cr
#' \code{30} = rotated BB8 copula (90 degrees) \cr
#' \code{33} = rotated Clayton copula (270 degrees) \cr
#' \code{34} = rotated Gumbel copula (270 degrees) \cr
#' \code{36} = rotated Joe copula (270 degrees) \cr
#' \code{37} = rotated BB1 copula (270 degrees) \cr
#' \code{38} = rotated BB6 copula (270 degrees) \cr
#' \code{39} = rotated BB7 copula (270 degrees) \cr
#' \code{40} = rotated BB8 copula (270 degrees) \cr
#' \code{104} = Tawn type 1 copula \cr
#' \code{114} = rotated Tawn type 1 copula (180 degrees) \cr
#' \code{124} = rotated Tawn type 1 copula (90 degrees) \cr
#' \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
#' \code{204} = Tawn type 2 copula \cr
#' \code{214} = rotated Tawn type 2 copula (180 degrees) \cr
#' \code{224} = rotated Tawn type 2 copula (90 degrees) \cr
#' \code{234} = rotated Tawn type 2 copula (270 degrees) \cr
#' @param param1 the first parameter for the bivariate copula
#' @param param2 the second parameter for the bivariate copula; if there is only one parameter, set param2 to 0 
#' 
#' @return log-likelihood vector for each pair of data
#' 
#' @examples
#' set.seed(100)
#' rho = .8
#' N = 100
#' censoring_rate = 0.5
#' u = VineCopula::BiCopSim(N, 1, rho, par2 = 0, obj = NULL, check.pars = TRUE)
#' x2 = round(u[, 2])
#' cdf = ecdf(x2)(x2)
#' cdf_unique = c(0, sort(unique(cdf)))
#' cdf_minus = cdf_unique[match(cdf, cdf_unique) - 1]
#' mixed_copula_llk(u[, 1], u[, 1], cdf, cdf_minus, FALSE, TRUE, 1, rho)
#'
#' @import VineCopula
#' @export

mixedcopulallk = function(u_plus, u_minus, v_plus, v_minus, is_disc_u, is_disc_v, fam, param1, param2=0) {
	u_plus[u_plus >= .999] = .999
	u_minus[u_minus >= .999] = .998
	v_plus[v_plus >= .999] = .999
	v_minus[v_minus >= .999] = .998
	u_plus[u_plus <= .001] = .002
	u_minus[u_minus <= .001] = .001
	v_plus[v_plus <= .001] = .002
	v_minus[v_minus <= .001] = .001
	
	if (is_disc_u == FALSE & is_disc_v == FALSE) {
		copula_llk = log(BiCopPDF(u_plus, v_plus, fam, param1, param2))
	} else if (is_disc_u == TRUE & is_disc_v == FALSE) {
		# condition on second variable (not discrete)
		numerator = BiCopHfunc2(u_plus, v_plus, fam, param1, param2) - BiCopHfunc2(u_minus, v_plus, fam, param1, param2)
		copula_llk = log(numerator / (u_plus - u_minus))
	} else if (is_disc_u == FALSE & is_disc_v == TRUE) {
		# condition on first variable (not discrete)
		numerator = BiCopHfunc2(v_plus, u_plus, fam, param1, param2) - BiCopHfunc2(v_minus, u_plus, fam, param1, param2)
		copula_llk = log(numerator / (v_plus - v_minus))
	} else if (is_disc_u == TRUE & is_disc_v == TRUE) {
		numerator = BiCopCDF(u_plus, v_plus, fam, param1, param2) - BiCopCDF(u_minus, v_plus, fam, param1, param2) - BiCopCDF(u_plus, v_minus, fam, param1, param2) + BiCopCDF(u_minus, v_minus, fam, param1, param2)
		copula_llk = log(numerator / (u_plus - u_minus) / (v_plus - v_minus))
	}
	return(copula_llk)
}