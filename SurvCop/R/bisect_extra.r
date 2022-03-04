#' Bisection method for finding a unique root
#' 
#' @description
#' Bisection method for a function which is monotone increasing in \eqn{[0,1]}.
#' Find root fun(v)=0 where fun is a function that returns a vector.
#' 
#' @param fun vector of strictly increasing functions on [0,1] which cross 0 
#' @param x0 vector of values where fun<0
#' @param x1 vector of values where fun>0
#' @param tol tolerance on consecutive iterates for stopping, default is 1.e-05
#' @param niter maximum number of iterations, default is 18
#' @param debug debug flag for intermediate steps, default is FALSE
#'
#' @return vector of roots 
#'
#' @examples
#' out <- bisect(function(x) sin(x)-c(.25,.5), c(0,0), c(1,1), niter=20, debug=TRUE)

bisect_extra=function(fun, x0, x1, tol=1e-05, niter=18, debug=FALSE, plus, minus) { 
	v0=x0; v1=x1
	diff=1; 
	f0=fun(v0,plus, minus); 
	f1=fun(v1,plus, minus)
	iter=0
	while(diff>tol & iter<=niter)
	{ v2=(v0+v1)/2
		f2=fun(v2,plus, minus)
		ii=(f2<0) 
		v0[ii]=v2[ii]; v1[!ii]=v2[!ii]
		diff=max(v1-v0)
		iter=iter+1
		if(debug) cat(iter,v2,"\n")
	}
	return(v2)
}