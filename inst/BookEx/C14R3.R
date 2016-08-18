## Initialization
MCMCsize <- 1e3
BurnIn <- trunc(MCMCsize / 2)
omega3ID <- omega3RW1 <- omega3RW2 <- rep(NA, MCMCsize)
omega3ID <- omega3RW1 <- omega3RW2 <- 0.5
acID <- acRW1 <- acRW2 <- 0
## MH: independence sampler
for(i in 2:MCMCsize){
    y <- rtrunc(1, spec = "exp",
                a = 0.0, b = 1.0, rate = 2)
    alpha <- min(1, (U1DU(y) / U1DU(omega3ID[i - 1])) *
                     (dtrunc(omega3ID[i - 1], spec = "exp",
                             a = 0.0, b = 1.0, rate = 2) /
                      dtrunc(y, spec = "exp",
                             a = 0.0, b = 1.0, rate = 2)))
    u <- runif(1)
    if(u < alpha){
        omega3ID[i] <- y
        acID <- acID + 1
    } else {
        omega3ID[i] <- omega3ID[i - 1]
    }
}
acrID <- acID / MCMCsize * 100
ar1ID <- ar(omega3ID, order.max = 1)$ar
omega3IDHat <- mean(omega3ID[-c(1:BurnIn)]) * 100
## MH: random walk sampler 
for(i in 2:MCMCsize){
    y1 <- omega3RW1[i - 1] + 0.5 * rnorm(1)
    if(y1 > 1 || y1 < 0) y1 <- -Inf
    y2 <- omega3RW2[i - 1] + 0.01 * rnorm(1)
    if(y2 > 1 || y2 < 0) y2 <- -Inf
    alpha1 <- min(1, U1DU(y1) / U1DU(omega3RW1[i - 1]))
    alpha2 <- min(1, U1DU(y2) / U1DU(omega3RW2[i - 1]))
    u <- runif(1)
    if(u < alpha1){
        omega3RW1[i] <- y1
        acRW1 <- acRW1 + 1
    } else {
        omega3RW1[i] <- omega3RW1[i - 1]
    }
    if(u < alpha2){
        omega3RW2[i] <- y2
        acRW2 <- acRW2 + 1
    } else {
        omega3RW2[i] <- omega3RW2[i - 1]
    }
}
acrRW1 <- acRW1 / MCMCsize * 100
ar1RW1 <- ar(omega3RW1, order.max = 1)$ar
omega3RW1Hat <- mean(omega3RW1[-c(1:BurnIn)]) * 100
acrRW2 <- acRW2 / MCMCsize * 100
ar1RW2 <- ar(omega3RW2, order.max = 1)$ar
omega3RW2Hat <- mean(omega3RW2[-c(1:BurnIn)]) * 100
