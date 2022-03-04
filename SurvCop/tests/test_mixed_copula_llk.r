library(SurvCop)

set.seed(100)
rho = .8
N = 100
censoring_rate = 0.5
u = VineCopula::BiCopSim(N, 1, rho, par2 = 0, obj = NULL, check.pars = TRUE)
x2 = round(u[, 2])
cdf = ecdf(x2)(x2)
cdf_unique = c(0, sort(unique(cdf)))
cdf_minus = cdf_unique[match(cdf, cdf_unique) - 1]

print(mixed_copula_llk(u[, 1], u[, 1], cdf, cdf_minus, FALSE, TRUE, 1, rho))