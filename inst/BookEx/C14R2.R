MCsim <- function(x0 = 0, n, theta0, theta1){
    ans <- vector()
    length(ans) <- n + 1
    ans[1] <- x0
    for(i in 2:(n + 1)){
        ans[i] <- theta0 + theta1 * ans[i-1] + rnorm(1)
        }
    ans
}
## Parameter settings
theta0 <- 2
theta1 <- 0.5
N <- 10000
x01 <- 14
x02 <- -10
## Markov Chains
mc1 <- MCsim(x0 = x01, n = N, theta0 = theta0, theta1 = theta1)
mc2 <- MCsim(x0 = x02, n = N, theta0 = theta0, theta1 = theta1)
## Progression of Markov Chains
plot(mc1[1:100], type = "l", ylim = range(cbind(mc1, mc2)),
     xlab = "", ylab = "X", main = "Progression of first-order Markov Chain")
lines(mc2[1:100], col = "red")
## Expected value of stationarity distribution and estimates
EfPop <- theta0 / (1 - theta1) 
m <- trunc(N / 2) ## burn-in
EfEst1 <- mean(mc1[-c(1:m)])
c(EfPop, EfEst1)
## Standard error of estimate for first MC
ar1Est <- ar(mc1, order.max = 1)$ar
se1 <- sqrt(var(mc1) / N * (1 + ar1Est) / (1 - ar1Est))
c(EfEst1 - 2 * se1, EfEst1, EfEst1 + 2 * se1)
## Variance of stationarity distribution and estimate
VarfPop <- 1 / (1 - theta1^2)
VarfEst1 <- 1 / (1 - ar1Est^2)
c(VarfPop, VarfEst1)
