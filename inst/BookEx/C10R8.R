library(cccp)
## Setting of parameters
Nassets <- ncol(rzoo)
Nobs <- nrow(rzoo)
mu <- colMeans(rzoo)
S <- cov(rzoo)
SR <- sqrm(S)
delta <- sqrt(qchisq(0.9, Nassets))
## Determining feasible risk aversion
SigMax <- max(colSds(rzoo))
SigMin <- min(colSds(rzoo))
ra <- seq(SigMin * 1.1, SigMax * 0.9, length.out = 10) / SigMax
## Initializing objects for MV and robust counterpart results
RCans <- MVans <- matrix(NA,
                         nrow = 10,
                         ncol = Nassets + 2)
## Computing points on efficient frontier and allocations
for(i in 1:10){
    ## minimum-variance
    wmv <- PMV(SRoot = SR, mu = mu, SigTerm = SigMin / ra[i])
    MVans[i, ] <- c(sqrt(t(wmv) %*% S %*% wmv),
                     crossprod(mu, wmv),
                     wmv)
    ## robust counterpart
    theta <- ra[i] + (1 - ra[i]) * delta / sqrt(Nobs)
    wrc <- PMV(SRoot = SR, mu = mu, SigTerm = SigMin / theta)
    RCans[i, ] <- c(sqrt(t(wrc) %*% S %*% wrc),
                     crossprod(mu, wrc),
                     wrc)
}
