#' van der Waerden correlation when first variable is continuous and second variable is right-censored
#' 
#' @description
#' van der Waerden correlation or correlation of normal scores when first variable is continuous and second variable is right-censored
#' 
#' @param dat1 variable 1 that is continuous, vector of length n
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
#' vWcor(x[,1], x[,2], censor_status)
#'
#' @import survival
#' @export

vWcor_continuous = function(dat1, dat2, censor_status, iprint=0, iterlim=100) {
	n = length(dat1)
	# convert variable 1 to N(0,1) via rank transform
	rank1 = rank(dat1)
	z1 = qnorm((rank1 - .5) / n)
	
	rank2 = numeric(n)
	# Kaplan-Meier empirical survival function
	# use (cdf[x+]+cdf[x-])/2 for all cases, censored or not
	F_2 = survival::survfit(Surv(dat2, 1 - censor_status) ~ 1)

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
		if( abs(rho)>=0.999) return(1.e10)
		nllk = 0
		if(ncens<n) { 
			tem=dbvn(cbind(z1[!icens],z2[!icens]),rho,ilog=T)
			tem[is.na(tem)]=-30
			nllk = nllk-sum(tem)
		}
		if(ncens>0) { 
			tem2=1-pnorm((z2[icens]-rho*z1[icens])/sqrt(1-rho^2))
			tem2[tem2<=0]=1.e-10
			nllk = nllk- sum(log(tem2)) 
		}
		return(nllk)
	}

	# starting value is spearman rho
	start= cor(dat1,dat2,method="spearman")
	rho_censored = nlm(f = nllk_censored, p = start, print.level = iprint, iterlim=iterlim)$estimate
	return(rho_censored)
}