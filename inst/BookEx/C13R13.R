## Evaluation of backtest
W <- list(WEP, WMD, WND)
FxRetBack <- FxRet[Eperiods, ]
Equity <- function(x, Eperiods, FxRet){
    WTs <- timeSeries(x, charvec = Eperiods)
    WTsL1 <- lag(WTs, k = 1)
    RetFacPR <- rowSums(FxRetBack / 100 * WTsL1) + 1
    RetFacPR[1] <- 100
    timeSeries(cumprod(RetFacPR), Eperiods)
}
WealthBack <- lapply(W, Equity, Eperiods = Eperiods, FxRet = FxRet)
## plot of wealth trajectories
ypretty <- unlist(lapply(WealthBack, function(x) pretty(range(x))))
ylims <- c(min(ypretty), max(ypretty))
plot(WealthBack[[1]], lwd = 2,
     xlab = "", ylab = "Index",
     main = "",
     ylim = ylims)
lines(WealthBack[[2]], lty = 2, lwd = 2, col = "blue")
lines(WealthBack[[3]], lty = 3, lwd = 2, col = "red")
legend("topleft", legend = c("EP", "Market", "Normal"), lty = 1:3,
       lwd = 2, col = c("black", "blue", "red"))
## portfolio performance/risk measures
PerfM <- function(x){
    EPRet <- returns(x, method = "discrete", percentage = FALSE,
                     trim = TRUE)
    EPRA <- Return.annualized(EPRet[, 1], scale = 52) * 100
    EPSD <- StdDev.annualized(EPRet[, 1], scale = 52) * 100
    EPSRA <- SharpeRatio.annualized(EPRet[, 1], scale = 52)
    EPES <- ES(EPRet[, 1], p = 0.95,
               method = "modified", clean = "boudt") * -100
    EPMDD <- maxDrawdown(EPRet[, 1]) * 100
    EPPerfW <- rbind(EPRA, EPSD, EPSRA, EPES, EPMDD)
    rownames(EPPerfW) <- c("Return (annual)", "Risk (annual, SD)",
                           "Sharpe Ratio", "CVaR (modified, 95%)", "Max Draw Down")
    EPPerfW
}
PerfBack <- lapply(WealthBack, PerfM)
PerfM <- matrix(unlist(PerfBack), ncol = 3)
rownames(PerfM) <- c("Return (annual)", "Risk (annual, SD)",
                     "Sharpe Ratio", "CVaR (modified, 95%)", "Maximum Draw Down")
colnames(PerfM) <- c("EP", "Market", "Normal")
PerfM
