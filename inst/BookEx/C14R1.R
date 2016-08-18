## Utility function
U1 <- function(x, mu, risk, lambda){
  lambda * x * mu - (1 - lambda) * risk * x^2
}
## Unscaled density
U1DU <- function(x, mu = 5, risk = 16, lambda = 0.5, nu = 1){
  exp(nu * U1(x = x, mu = mu, risk = risk, lambda = lambda))
}
## Constant of density
Z <- integrate(U1DU, lower = 0, upper = 1)$value
## Classical MC integration for expected value of scaled U1DU
MCsize <- 1e3
u1 <- runif(MCsize)
omega1 <- u1 * U1DU(u1) / Z * 100
omega1Hat <- mean(omega1)
## Graphical presentation of recursive means with bands
omega1RecMean <- cumsum(omega1) / (1:MCsize)
omega1RecErr <- sqrt(cumsum((omega1  - omega1RecMean)^2)) / (1:MCsize)
plot(omega1RecMean, type = "l",
     ylim = c(0, 100), ylab = expression(omega), xlab = "",
     main = "Recursive estimates for allocation to risky asset")
lines(omega1RecMean + 2 * omega1RecErr, col = "red", lty = 2)
lines(omega1RecMean - 2 * omega1RecErr, col = "red", lty = 2)
## Random draws from U1DU by accept-reject method
## Candidate distribution is a truncated Exponential with lambda  = 2 and k = 1.5
library(truncdist)
k <- 1.5
y2 <- rtrunc(MCsize, spec = "exp",
             a = 0.0, b = 1.0, rate = 2)
u2 <- runif(MCsize)
AcceptRejectIdx <- which(u2 <= (U1DU(y2) /
                                (k * dtrunc(y2, spec = "exp", a = 0, b = 1, rate = 2))))
omega2 <- y2[AcceptRejectIdx]
omega2Hat <- mean(omega2) * 100
