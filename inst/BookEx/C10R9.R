## Equivalent weighting parameter for theta
theta2lambda <- function(theta, Nassets, Nobs, level){
  delta <- sqrt(qchisq(level, Nassets))
  lambda <- theta / (1 + theta * (delta / sqrt(Nobs)))
  lambda
}
## robust allocation and equivalent
## mean-variance allocation
theta <- 0.7
wrc <- PMV(SRoot = SR, mu = mu, SigTerm = SigMin / theta)
## RC point on efficient frontier
rceq <- c(sqrt(t(wrc) %*% S %*% wrc), crossprod(mu, wrc))
## Equivalent risk weighting
rweq <- theta2lambda(theta, Nassets = Nassets, Nobs = Nobs, level = 0.9)
## Equivalent MV point on efficient frontier
wmv <- PMV(SRoot = SR, mu = mu, SigTerm = SigMin / rweq)
mveq <- c(sqrt(t(wmv) %*% S %*% wmv), crossprod(mu, wmv))
