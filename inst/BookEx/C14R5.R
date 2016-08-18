## Loading of packages
library(FRAPO)
library(adaptMCMC)
library(MCMCpack)
library(mcmc)
library(rstan)
library(Rsolnp)
## Computation of equity returns
data(MultiAsset)
Assets <- timeSeries(MultiAsset, charvec = rownames(MultiAsset))
R <- 100 * ((1 + returns(Assets, method = "discrete", percentage = FALSE))^12 - 1)
#R <- returns(Assets, method = "discrete", percentage = FALSE)
Requity <- R[, c("GSPC", "RUA", "GDAXI", "FTSE")]
## Parameter settings/initialization
muSample <- colMeans(Requity)
sigmaSample <- cov(Requity) 
(targetR <- mean(muSample))
targetR <- 10
riskA <- 3
Nt <- 60
convR <- sqrt(nrow(Requity))
N <- ncol(Requity)
##
## Maximizing expected utility
##
## objective function to be minimized
f0 <- function(pars, mu, Sigma, Lambda, L, iperiod){
    uval <- U(pars, mu, Sigma, Lambda, L, iperiod)
    -1.0 * uval
}
## budget constraint for optimization
budget <- function(pars, mu, Sigma, Lambda, L, iperiod){
    sum(pars)
}
## initial point
w0 <- rep(1/N, N)
## Computing MEU allocation
optMEU <- solnp(pars = w0, fun = f0, eqfun = budget, eqB = 1.0,
                LB = rep(0, N), mu = muSample, Sigma = sigmaSample,
                Lambda = riskA, L = targetR, iperiod = Nt)
wMEU <- optMEU$pars * 100
UMEU <- U(wMEU / 100, mu = muSample, Sigma = sigmaSample,
          Lambda = riskA, L = targetR, iperiod = Nt)
