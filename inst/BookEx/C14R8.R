## Analyzing results
## Computing distances
MUD <- PUD <- matrix(NA, nrow = Draws, ncol = SS)
for(i in 1:SS){
    MUD[, i] <- (UMEU - MU[, i]) / abs(UMEU)
    PUD[, i] <- (UMEU - PU[, i]) / abs(UMEU)
}
## Graph
PuDMean <- apply(PUD, 2, function(x) mean(x))
MuDMean <- apply(MUD, 2, function(x) mean(x))
PuDMin <- apply(PUD, 2, function(x) min(x))
MuDMin <- apply(MUD, 2, function(x) min(x))
PuDMax <- apply(PUD, 2, function(x) max(x))
MuDMax <- apply(MUD, 2, function(x) max(x))
ylims <- range(na.omit(c(PuDMax, MuDMax, PuDMin, MuDMin)))
plot(cbind(1:SS, PuDMean), type = "p", col = "blue", pch = 17, cex = 1.2,
     ylim = ylims, ylab = "Relative deviations from 'true' utility",
     xlab = "Sample Sizes", axes = FALSE) 
points(1:SS, PuDMax, type = "p", col = "blue", cex = 1.2, pch = 22)
points(1:SS, PuDMin, type = "p", col = "blue", cex = 1.2, pch = 22)
points(1:SS, MuDMean, type = "p", col = "darkred", cex = 1.2, pch = 19)
points(1:SS, MuDMax, type = "p", col = "darkred", cex = 1.2, pch = 22)
points(1:SS, MuDMin, type = "p", col = "darkred", cex = 1.2, pch = 22)
arrows(x0 = 1:SS, y0 = PuDMin, y1 = PuDMax, code = 3, col = "blue",
       length = 0.0, angle = 90, lwd = 2)
arrows(x0 = 1:SS, y0 = MuDMin, y1 = MuDMax, code = 3, col = "darkred",
       length = 0.0, angle = 90, lwd = 1.0)
axis(1, at = 1:SS, Samples)
axis(2, pretty(ylims), pretty(ylims), las = 2) 
box()
legend("topright", legend = c("PU mean deviation",
                       "PU min/max deviation",
                       "MU mean deviation",
                       "MU min/max deviation"),
       pch = c(17, 22, 19, 22),
       col = c(rep("blue", 2), rep("darkred", 2)))
abline(h = 0, col = "gray")
