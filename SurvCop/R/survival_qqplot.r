#' generate the Q-Q plot of the empirical quantiles values versus the quantiles from a parametric model for a right-censored variable
#' 
#' @description
#' generate the Q-Q plot of the empirical quantiles values versus the quantiles from a parametric model for a right-censored variable
#' 
#' @param y n by 1 numeric vector of the values of the right-censored variable
#' @param censor_status n by 1 logical vector, ith position TRUE if udata[i,d] is right-censored
#' @param qfunc function to get the quantile from a parametric model, parameters of the parametric model are provided as extra arguments at the end of the function call
#
#' @return the Q-Q plot of the empirical quantiles values versus the quantiles from a parametric model for a right-censored variable
#'
#' @import survival
#' @export 

survival_qqplot = function(y, censor_status, qfunc, ...) {
	km_fit = survfit(Surv(y, 1 - censor_status) ~ 1)
	emp_cdf = 1 - km_fit$surv
	emp_cdf = (emp_cdf + c(0, emp_cdf[-1])) / 2
	y_p = qfunc(emp_cdf, ...)
	plot(km_fit$time~y_p, ylab = 'empirical quantiles', xlab = 'parametric quantiles', main = 'QQ plot')
	abline(c(0, 1))
}
