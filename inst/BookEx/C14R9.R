## Mean allocations of PU and MU
PUWmean <- MUWmean <- matrix(NA, nrow = SS, ncol = N)
for(i in 1:SS){
    PUWmean[i, ] <- colMeans(PUW[, , i])
    MUWmean[i, ] <- colMeans(MUW[, , i])
}
## Min of PU and MU allocations
PUWmin <- MUWmin <- matrix(NA, nrow = SS, ncol = N)
for(i in 1:SS){
    PUWmin[i, ] <- apply(PUW[, , i], 2, min)
    MUWmin[i, ] <- apply(MUW[, , i], 2, min)
}
## Max of PU and MU allocations
PUWmax <- MUWmax <- matrix(NA, nrow = SS, ncol = N)
for(i in 1:SS){
    PUWmax[i, ] <- apply(PUW[, , i], 2, max)
    MUWmax[i, ] <- apply(MUW[, , i], 2, max)
}
## Range of PU and MU allocations
PUWrange <- PUWmax - PUWmin
MUWrange <- MUWmax - MUWmin
rownames(PUWmean) <- paste("Sample", Samples, sep = "-")
colnames(PUWmean) <- colnames(Requity)
rownames(MUWmean) <- rownames(PUWmean)
colnames(MUWmean) <- colnames(Requity)
rownames(PUWrange) <- paste("Sample", Samples, sep = "-")
colnames(PUWrange) <- colnames(Requity)
rownames(MUWrange) <- rownames(PUWrange)
colnames(MUWrange) <- colnames(Requity)
