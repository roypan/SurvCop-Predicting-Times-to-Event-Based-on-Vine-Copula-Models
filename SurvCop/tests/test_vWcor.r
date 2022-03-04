library(SurvCop)

set.seed(100)
rho = .8
N = 1000
censoring_rate = 0.5
u = VineCopula::BiCopSim(N, 1, rho, par2 = 0, obj = NULL, check.pars = TRUE)
auxiliary = pnorm(rnorm(N, -sqrt(2) * qnorm(censoring_rate)))
censor_status = auxiliary < u[, 2]
u[as.logical(censor_status), 2] = auxiliary[as.logical(censor_status)]
x = qnorm(u)
print(vWcor_continuous(x[,1], x[,2], censor_status))

set.seed(100)
rho = .8
N = 1000
censoring_rate = 0.5
u = VineCopula::BiCopSim(N, 1, rho, par2 = 0, obj = NULL, check.pars = TRUE)
auxiliary = pnorm(rnorm(N, -sqrt(2) * qnorm(censoring_rate)))
censor_status = auxiliary < u[, 2]
u[as.logical(censor_status), 2] = auxiliary[as.logical(censor_status)]
x = qnorm(u)
x[, 1] = round(x[, 1])
print(vWcor_discrete(x[,1], x[,2], censor_status))