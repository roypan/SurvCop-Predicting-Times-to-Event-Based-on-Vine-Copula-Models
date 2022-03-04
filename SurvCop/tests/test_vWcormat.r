library(SurvCop)

set.seed(100)
d = 5
censoring_rate = 0.5
N = 1000

# define 5-dimensional R-vine tree structure matrix
Matrix = c(5, 1, 2, 3, 4,
		   0, 4, 1, 2, 3,
		   0, 0, 2, 1, 3,
		   0, 0, 0, 3, 1,
		   0, 0, 0, 0, 1)
Matrix = matrix(Matrix, d, d)

# define R-vine pair-copula family matrix
family = c(0, 1, 1, 1, 1,
		   0, 0, 1, 1, 1,
		   0, 0, 0, 1, 1,
		   0, 0, 0, 0, 1,
		   0, 0, 0, 0, 0)
family = matrix(family, d, d)

par = matrix(0.6, nrow = d, ncol = d)
par2 = matrix(0, nrow = d, ncol = d)

RVM = RVineMatrix(Matrix = Matrix, family = family, par = par, par2 = par2)

u = RVineSim(N, RVM)
auxiliary = pnorm(rnorm(N, -sqrt(2) * qnorm(censoring_rate)))
censor_status = auxiliary < u[, d]
u[as.logical(censor_status), d] = auxiliary[as.logical(censor_status)]
is_disc = c(FALSE, TRUE, TRUE, FALSE, FALSE)
u[, is_disc] = round(u[, is_disc])

print(vWcormat(u, is_disc, censor_status))
