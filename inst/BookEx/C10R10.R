## Efficient Frontier
## determining ranges
xlims <- c(SigMin * 0.9,
           max(cbind(MVans[, 1], RCans[, 1])))
ylims <- c(min(cbind(MVans[, 2], RCans[, 2])),
           max(cbind(MVans[, 2], RCans[, 2])))
## plotting efficient frontier for MV
plot(MVans[, 1], MVans[, 2], type = "l",
     xlim = xlims,
     ylim = ylims,
     xlab = expression(sigma),
     ylab = expression(mu))
## Superimposing points
for(i in 1:nrow(MVans)){
    points(x = MVans[i, 1], y = MVans[i, 2], col = "blue", pch = 19)
    points(x = RCans[i, 1], y = RCans[i, 2], col = "red", pch = 19)
}
## Superimposing equivalence points
points(x = rceq[1], y = rceq[2], col = "darkgreen", bg = "darkgreen", pch = 23)
points(x = mveq[1], y = mveq[2], col = "green", bg = "green", pch = 23)
## Legend
legend("topleft", legend = c("Efficient Frontier", "MV points", "RC points",
                      expression(paste("RC allocation with ", theta == 0.7)),
                      "Equivalent MV-Portfolio"),
       lty = c(1, NA, NA, NA, NA), pch = c(NA, 19, 19, 23, 23),
       pt.bg = c(NA, NA, NA, "darkgreen", "orange"),
       col = c("black", "blue", "red", "darkgreen", "orange"))
