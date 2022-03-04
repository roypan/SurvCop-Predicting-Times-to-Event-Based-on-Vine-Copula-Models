library(SurvCop)

SurvCop:::BivCopPDF(.5,.3,2,.6,6.6)
BiCopPDF(.5,.3,2,.6,6.6)
SurvCop:::BivCopCCDF(.5,.3,2,.6,6.6)
BiCopHfunc2(.5,.3,2,.6,6.6)

SurvCop:::mixedcopulallk(.5,.4,.3,.2, F, F, 7, 2, 5)
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, F, F, 7, 2, 5)
SurvCop:::mixedcopulallk(.5,.4,.3,.2, F, T, 7, 2, 5)
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, F, T, 7, 2, 5)
SurvCop:::mixedcopulallk(.5,.4,.3,.2, T, F, 7, 2, 5)
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, T, F, 7, 2, 5)
SurvCop:::mixedcopulallk(.5,.4,.3,.2, T, T, 7, 2, 5)
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, T, T, 7, 2, 5)

SurvCop:::mixedcopulallk(.5,.4,.3,.2, F, F, 17, 2, 5)
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, F, F, 17, 2, 5)
SurvCop:::mixedcopulallk(.5,.4,.3,.2, F, T, 17, 2, 5)
#[1] NaN
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, F, T, 17, 2, 5)
SurvCop:::mixedcopulallk(.5,.4,.3,.2, T, F, 17, 2, 5)
#[1] NaN
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, T, F, 17, 2, 5)
SurvCop:::mixedcopulallk(.5,.4,.3,.2, T, T, 17, 2, 5)
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, T, T, 17, 2, 5)

SurvCop:::mixedcopulallk(.5,.4,.3,.2, F, F, 5, 6)
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, F, F, 5, 6)
SurvCop:::mixedcopulallk(.5,.4,.3,.2, F, T, 5, 6)
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, F, T, 5, 6)
SurvCop:::mixedcopulallk(.5,.4,.3,.2, T, F, 5, 6)
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, T, F, 5, 6)
SurvCop:::mixedcopulallk(.5,.4,.3,.2, T, T, 5, 6)
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, T, T, 5, 6)

SurvCop:::mixedcopulallk(.5,.4,.3,.2, F, F, 3, .5)
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, F, F, 3, .5)
SurvCop:::mixedcopulallk(.5,.4,.3,.2, F, T, 3, .5)
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, F, T, 3, .5)
SurvCop:::mixedcopulallk(.5,.4,.3,.2, T, F, 3, .5)
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, T, F, 3, .5)
SurvCop:::mixedcopulallk(.5,.4,.3,.2, T, T, 3, .5)
SurvCop:::mixed_copula_llk(.5,.4,.3,.2, T, T, 3, .5)

d = 5
n=100
Matrix = c(5, 1, 2, 3, 4,
		   0, 4, 1, 2, 3,
		   0, 0, 2, 1, 3,
		   0, 0, 0, 3, 1,
		   0, 0, 0, 0, 1)
Matrix = matrix(Matrix, d, d)

family = c(0, 5, 5, 5, 5,
		   0, 0, 5, 1, 4,
		   0, 0, 0, 1, 4,
		   0, 0, 0, 0, 4,
		   0, 0, 0, 0, 0)
family = matrix(family, d, d)

par1 = c(0, 4, 5, 3, 4,
		 0, 0, 5, .5, 4,
		 0, 0, 0, .4, 3,
		 0, 0, 0, 0, 5,
		 0, 0, 0, 0, 0)
par1 = matrix(par1, d, d)

set.seed(10)
data_cdf_plus = matrix(runif(n*d, .4, .95), n, d)
data_cdf_minus = data_cdf_plus - .1
censor_status = as.logical(rbinom(n,1,.5))
is_disc = c(FALSE, TRUE, TRUE, FALSE, FALSE)

SurvCop:::survival_nllk(paramvec = par1[d:2, 1], param1_complete=par1[d:2,d:2], param2_complete=matrix(0,d-1,d-1), fam=family[d:1,d:1], A=Matrix[d:1,d:1], data_cdf_plus, data_cdf_minus, censor_status, is_disc)

SurvCop:::survivalnllk(paramvec = par1[d:2, 1], param1_complete=par1[d:2,d:2], param2_complete=matrix(0,d-1,d-1), fam=family[d:1,d:1], A=Matrix[d:1,d:1], data_cdf_plus, data_cdf_minus, censor_status, is_disc)

d = 5
n=100
Matrix = c(5, 1, 2, 3, 4,
		   0, 4, 1, 2, 3,
		   0, 0, 2, 1, 3,
		   0, 0, 0, 3, 1,
		   0, 0, 0, 0, 1)
Matrix = matrix(Matrix, d, d)

family = c(0, 16, 7, 3, 4,
		   0, 0, 17, 6, 4,
		   0, 0, 0, 13, 14,
		   0, 0, 0, 0, 14,
		   0, 0, 0, 0, 0)
family = matrix(family, d, d)

par1 = c(0, 3, 3, 3, 8,
		 0, 0, 3, 3, 8,
		 0, 0, 0, 3, 9,
		 0, 0, 0, 0, 8,
		 0, 0, 0, 0, 0)
par1 = matrix(par1, d, d)

par2 = c(0, 3, 3, 3, 4,
		 0, 0, 3, 3, 4,
		 0, 0, 0, 3, 3,
		 0, 0, 0, 0, 5,
		 0, 0, 0, 0, 0)
par2 = matrix(par2, d, d)

set.seed(100)
data_cdf_plus = matrix(runif(n*d, .3, .35), n, d)
data_cdf_minus = data_cdf_plus - .1
censor_status = as.logical(rbinom(n,1,.5))
is_disc = c(TRUE, TRUE, TRUE, TRUE, FALSE)

SurvCop:::survival_nllk(paramvec = par1[d:2, 1], param1_complete=par1[d:2,d:2], param2_complete=par2[d:2,d:2], fam=family[d:1,d:1], A=Matrix[d:1,d:1], data_cdf_plus, data_cdf_minus, censor_status, is_disc)

SurvCop:::survivalnllk(paramvec = par1[d:2, 1], param1_complete=par1[d:2,d:2], param2_complete=par2[d:2,d:2], fam=family[d:1,d:1], A=Matrix[d:1,d:1], data_cdf_plus, data_cdf_minus, censor_status, is_disc)

d = 5
n=1000
Matrix = c(5, 1, 2, 3, 4,
		   0, 4, 1, 2, 3,
		   0, 0, 2, 1, 3,
		   0, 0, 0, 3, 1,
		   0, 0, 0, 0, 1)
Matrix = matrix(Matrix, d, d)

family = c(0, 5, 5, 5, 5,
		   0, 0, 5, 1, 4,
		   0, 0, 0, 1, 4,
		   0, 0, 0, 0, 4,
		   0, 0, 0, 0, 0)
family = matrix(family, d, d)

par1 = c(0, 4, 5, 3, 4,
		 0, 0, 5, .5, 4,
		 0, 0, 0, .4, 3,
		 0, 0, 0, 0, 5,
		 0, 0, 0, 0, 0)
par1 = matrix(par1, d, d)

set.seed(10)
data_cdf_plus = matrix(runif(n*(d-1), .4, .95), n, d-1)
data_cdf_minus = data_cdf_plus - .1
is_disc = c(FALSE, TRUE, TRUE, FALSE, FALSE)
u = runif(n, .1, .9)

start_time <- Sys.time()
SurvCop:::survival_conditional_cdf(u, paramvec = par1[d:2, 1], param1_complete=par1[d:2,d:2], param2_complete=matrix(0,d-1,d-1), fam=family[d:1,d:1], A=Matrix[d:1,d:1], data_cdf_plus, data_cdf_minus, is_disc)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
SurvCop:::survivalconditionalcdf(u, paramvec = par1[d:2, 1], param1_complete=par1[d:2,d:2], param2_complete=matrix(0,d-1,d-1), fam=family[d:1,d:1], A=Matrix[d:1,d:1], data_cdf_plus, data_cdf_minus, is_disc)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
SurvCop:::survivalpredict(p=c(.25,.5,.75), paramvec = par1[d:2, 1], param1_complete=par1[d:2,d:2], param2_complete=matrix(0,d-1,d-1), fam=family[d:1,d:1], A=Matrix[d:1,d:1], data_cdf_plus, data_cdf_minus, is_disc, iprint=F)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
SurvCop:::survival_predict(p=c(.25,.5,.75), paramvec = par1[d:2, 1], param1_complete=par1[d:2,d:2], param2_complete=matrix(0,d-1,d-1), fam=family[d:1,d:1], A=Matrix[d:1,d:1], data_cdf_plus, data_cdf_minus, is_disc, iprint=F)
end_time <- Sys.time()
end_time - start_time

