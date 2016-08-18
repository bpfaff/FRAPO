## Backtest
for(i in 1:LengthBack){
  cat(paste("Backtestperiod:", Eperiods[i], "\n"))
  set.seed(i + 12345)
  fp <- window(FxRet, start = Speriods[i], end = Eperiods[i])
  ## Update Market Model
  MarketMod <- mstFit(series(fp))
  par <- MarketMod@fit$estimated
  M <- rmst(J, xi = c(par$beta), Omega = par$Omega,
            alpha = par$alpha, nu = par$nu)
  colnames(M) <- Anames
  mu2 <- colMeans(M)^2
  Vt <- t(M)
  Aeq <- rbind(Vt^2, rep(1, J))
  ## GARCH-model
  mfcst <- lapply(fp, CondVolaFcst)
  fcst <- matrix(unlist(mfcst), ncol = 1, byrow = TRUE)
  beq <- matrix(c(mu2 + fcst[, 1]^2, 1), ncol = 1)
  ## EP-optimization
  Ep <- ep(x0, Aeq, beq, pprior)
  ## EP for fixed tau = 0.5
  EpH <- 0.5 * Ep + 0.5 * pprior
  EPspec <- portfolioSpec()
  EPmom <- function(x, spec = NULL, ...){
      m <- cov.wt(x, drop(EpH))
      list("mu" = m$center, "Sigma" = m$cov)
  }
  setEstimator(EPspec) <- "EPmom"
  PSol <- tangencyPortfolio(data = as.timeSeries(M), spec = EPspec)
  WEP[i, ] <- getWeights(PSol)
  ## Portfolio based on market distribution
  PSolM <- tangencyPortfolio(data = as.timeSeries(M), spec = portfolioSpec())
  WMD[i, ] <- getWeights(PSolM)
  ## Portfolio based on normality assumption
  PSolN <- tangencyPortfolio(data = fp, spec = portfolioSpec())
  WND[i, ] <- getWeights(PSolN)
}
