#' van der Waerden correlation when first variable is discrete and second variable is right-censored
#' 
#' @description
#' van der Waerden correlation or correlation of normal scores when first variable is discrete and second variable is right-censored
#' 
#' @param dat1 variable 1 that is discrete, vector of length n
#' @param dat2 variable 2 that is right-censored, vector of length n
#' @param censor_status binary vector of length n, indicator of censored , =0 for uncensored, =1 for right-censored
#' @param iprint print_level for nlm minimization
#' @param iterlim number of iterations for nlm minimization
#' 
#' @return van der Waerden (normal score) correlation value in the interval (-1,1)
#' 
#' @examples
#' rho = .8
#' N = 1000
#' censoring_rate = 0.5
#' u = VineCopula::BiCopSim(N, 1, rho, par2 = 0, obj = NULL, check.pars = TRUE)
#' auxiliary = pnorm(rnorm(N, -sqrt(2) * qnorm(censoring_rate)))
#' censor_status = auxiliary < u[, 2]
#' u[as.logical(censor_status), 2] = auxiliary[as.logical(censor_status)]
#' x = qnorm(u)
#' x[, 1] = round(x[, 1])
#' vWcor(x[,1], x[,2], censor_status)
#' 
#' @import survival
#' @export

vWcor_discrete = function(dat1, dat2, censor_status, iprint=0, iterlim=100) {
	# bivariate normal negative log-likelihood when second variable is right-censored 
	# param is rho parameter in (-1,1) and a negative log-likelihood is returned
	n = length(dat1)
	# convert variable 1 to N(0,1) via rank transform
	F_1 = ecdf(dat1)(dat1)
	F_1_unique = sort(unique(F_1))
	F_1_unique = c(0, F_1_unique)
	z1_u = numeric(n)
	z1_l = numeric(n)
	for (i in 1:n) {
		ind = which(F_1_unique == F_1[i])
		z1_u[i] = qnorm(F_1_unique[ind])
		z1_l[i] = qnorm(F_1_unique[ind - 1])
	}
	# modification for Inf
	ulim=6
	z1_u[z1_u>= ulim]=qnorm(.999)
	z1_l[z1_l<= -ulim]=qnorm(.001)
	rank2 = numeric(n)
	# Kaplan-Meier empirical survival function
	# use (cdf[x+]+cdf[x-])/2 for all cases, censored or not
	F_2 = survfit(Surv(dat2, 1 - censor_status) ~ 1)
	for (i in 1:n) {
		ind = which.min(abs(dat2[i] - F_2$time))
		y=dat2[i]
		dat2_minus = dat2[dat2 < y] # find smaller times
		if (length(dat2_minus) == 0) { # if it is the smallest observation
			rank2[i] = (1-F_2$surv[ind])/2 # cdf[x-]=0, cdf[x]=1-surv[x]
		} else { # F_2$surv[ind2] > F_2$surv[ind]
			ind2 = which.min(abs(max(dat2_minus) - F_2$time))
			ind2 = ind2[1] # in case length is more than one 
			rank2[i] = (2 -F_2$surv[ind]-F_2$surv[ind2])/2 #(cdf[x]+cdf[x-])/2
		}
	}
	z2 = qnorm(rank2)
	icens = (censor_status == 1)
	ncens=sum(icens)

	# start of internal function for optimization	 
	nllk_censored = function(param) {
		rho = param
		sigma = matrix(c(1, rho, rho, 1), nrow = 2)
		if( abs(rho)>=0.999) return(1.e10)
		r1=sqrt(1-rho^2)
		nllk = 0
		# if dat2 is not censored
		if(ncens<n) {
			tem2=pnorm((z1_u[!icens]-rho*z2[!icens])/r1)
			tem1=pnorm((z1_l[!icens]-rho*z2[!icens])/r1)
			tem2[tem2>=1]=1-1.e-10
			tem1[tem1<=0]=1.e-10
			diff=tem2-tem1
			diff[diff<=0]=1.e10 # prevent NaN in intermediate steps of nlm
			nllk = nllk - sum(log(diff))
		}
		
		if(ncens>0) {
			for (i in 1:n) {
				if (icens[i]) {
					rect=pmvnorm(lower = c(z1_l[i], z2[i]), upper = c(z1_u[i], Inf), mean = c(0, 0), sigma = sigma)
					if(rect<=0) rect=1.e-10
					nllk = nllk - log(rect)
				}
			}
		}
		return(nllk)
	}
	
	# starting value is spearman rho
	start= cor(dat1,dat2,method="spearman")
	rho_censored = nlm(f=nllk_censored, p=start, print.level=iprint, iterlim=iterlim)$estimate
	return(rho_censored)
}
