##
## Auxilliary functions for utility
##
## complementary error function
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
## Downside
f1 <- function(x, Lambda, eta){
    (Lambda - 1) * exp(-x * eta^2) / sqrt(2 * pi) 
}
## Upside
f2 <- function(x, Lambda, eta){
    (1 + Lambda) / 2 - (Lambda - 1) / 2 * erfc(sqrt(x) * eta)
}
##
## Utility
##
U <- function(w, mu, Sigma, Lambda, L, iperiod){
    M <- drop(crossprod(w, mu) - L)
    S <- drop(sqrt(crossprod(w, Sigma) %*% w))
    eta <- M / (2 * S)
    ti <- 1:iperiod
    dt <- 1 / iperiod
    f1val <- f1(ti * dt, Lambda, eta)
    f2val <- f2(ti * dt, Lambda, eta)
    uval <- sum(dt * (ti * dt * M * f2val - sqrt(ti * dt) * S * f1val))
    uval
}
##
## Probabilistic re-interpretation of utility 
## (unscaled log-density)
##
PUL <- function(w, mu, Sigma, Lambda, L, iperiod, nu){
    if(any(w < 0)){
        return(-Inf)
    } else if(any(w >= 1)){
        return(-Inf)
    } else {
        w <- w / sum(w)
        uval <- U(w, mu, Sigma, Lambda, L, iperiod)
        return(nu * uval)
    }
}
