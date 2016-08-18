library(FRAPO)
library(fPortfolio)
## Loading of data set
data(MultiAsset)
## Creating timeSeries of prices and returns
pr <- timeSeries(MultiAsset, charvec = rownames(MultiAsset))
data <- returns(pr, methdo = "discrete", percentages = TRUE, trim = TRUE)
## Parameters / constant
NAssets <- ncol(pr)
ANames <- colnames(pr)
Sigma <- cov(data)
mu <- colMeans(data)
## Risk surface plot
hull <- markowitzHull(data, nFrontierPoints = 50)
grid <- feasibleGrid(hull, trace = FALSE)
divers <- bestDiversification(grid, trace = FALSE)
## Standard deviation of marginal risk contributions
mrc.sd <- function(data, weights){
    Sigma <- cov(data)
    a <- mrc(weights, Sigma)
    sd(a)
}
surf <- riskSurface(divers, FUN = "mrc.sd")
## Feasible portfolios with highest diversification ratio
allWeights <- attr(divers, "weights")
idx <- sort(unique(allWeights[, 1]))
dropt <- matrix(0, nrow = length(idx), ncol = 2)
idxRow <- 1:length(idx)
for(j in idx){
    w <- matrix(allWeights[allWeights[, 1] == j, -c(1, 2)], ncol = NAssets)
    divm <- vector()
    length(divm) <- nrow(w)
    for(i in 1:nrow(w)){
        divm[i] <- dr(w[i, ], Sigma)
    }
    divmidx <- which.max(divm)
    wopt <- w[divmidx, ]
    dropt[idxRow[j], ] <- c(crossprod(wopt, mu),
                            sqrt(crossprod(wopt, Sigma) %*% wopt))
}
