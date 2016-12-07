## Defining function for points on efficient frontier
PMV <- function(SRoot, mu, SigTerm,
                    optctrl = ctrl(trace = FALSE)){
  N <- nrow(SRoot)
  ## Portfolio risk constraint
  soc1 <- socc(F = SRoot, g = rep(0, N), d = rep(0, N), f = SigTerm)
  ## non-negativity constraint
  nno1 <- nnoc(G = -diag(N), h = rep(0, N))
  ## Budget constraint
  A1 <- matrix(rep(1, N), nrow = 1)
  b1 <- 1.0
  ## optimization
  ans <- cccp(q = -mu, A = A1, b = b1, cList = list(nno1, soc1),
              optctrl = optctrl)
  getx(ans)
}
