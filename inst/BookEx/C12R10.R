sdmrc <- surf$z
c104 <- which((sdmrc >= 10.35) & (sdmrc <= 10.45),
              arr.ind = TRUE)
w104 <- matrix(NA, nrow = nrow(c104), ncol = NAssets)
colnames(w104) <- ANames
for(i in 1:nrow(c104)){
    gidx <- which((allWeights[, 1] == c104[i, 1]) &
                  (allWeights[, 2] == c104[i, 2]), arr.ind = TRUE)
    w104[i, ] <- allWeights[gidx, -c(1, 2)]
}
## Computing standard deviations of mrc and standard deviation risk
sdmrc104 <- apply(w104, 1, function(x) sd(mrc(x, Sigma = Sigma)))
sdr104 <- apply(w104, 1, function(x) sqrt(crossprod(x, Sigma) %*% x)) * 100
## Grouping by asset class
wEquity <- w104[, 1:6]
wBonds <- w104[, 7:9]
wGold <- w104[, 10]
wEquity <- rowSums(wEquity)
wBonds <- rowSums(wBonds)
wAsset <- cbind(wEquity, wBonds, wGold) * 100
ans <- cbind(wAsset, sdmrc104, sdr104)
colnames(ans) <- c("Equity", "Bonds", "Gold", "StdDev. of MRC", "StdDev. Risk")
rownames(ans) <- 1:nrow(ans)
