#' generate the Q-Q plot of the empirical CDF values versus the CDF values from a parametric model for a right-censored variable
#' 
#' @description
#' generate the Q-Q plot of the empirical CDF values (on the y-axis) versus the CDF values from a parametric model (on the x-axis) for a right-censored variable
#' 
#' @param y n by 1 numeric vector of the values of the right-censored variable
#' @param censor_status n by 1 logical vector, ith position TRUE if udata[i,d] is right-censored
#' @param pfunc function to get the CDF from a parametric model, parameters of the parametric model are provided as extra arguments at the end of the function call
#
#' @return the Q-Q plot of the empirical CDF values versus the CDF values from a parametric model for a right-censored variable
#'
#' @import survival
#' @export 

survival_qqplot = function(y, censor_status, pfunc, ...) {
	km_fit = survfit(Surv(y, 1 - censor_status) ~ 1)
	emp_cdf = 1 - km_fit$surv
	uy_p = pfunc(km_fit$time, ...)
	plot(emp_cdf~uy_p, ylab = 'empirical CDF', xlab = 'parametric CDF', main = 'QQ plot')
	abline(c(0, 1))
}