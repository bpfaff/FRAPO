## Function for utility simulation
PuMuSim <- function(x){
    SmpL <- nrow(x)
    muS <- colMeans(x)
    sigmaS <- cov(x)
    PUans <- metrop(obj = PUL, initial = w0, nbatch = McmcLength, blen = 1,
                    scale = 1e-2, mu = muS, Sigma = sigmaS,
                    Lambda = riskA, L = targetR, iperiod = Nt, nu = sqrt(SmpL))
    PUW <- colMeans(PUans$batch)
    PUW <- PUW / sum(PUW)
    PU <- U(w = PUW, mu = muSample, Sigma = sigmaSample, Lambda = riskA, L = targetR,
                   iperiod = Nt)
    MUans <- solnp(pars = w0, fun = f0, eqfun = budget, eqB = 1.0, LB = rep(0, N),
                   mu = muS, Sigma = sigmaS, Lambda = riskA, L = targetR,
                   iperiod = Nt)
    MUW <- MUans$pars
    MU <- U(w = MUW, mu = muSample, Sigma = sigmaSample, Lambda = riskA, L = targetR,
                   iperiod = Nt)
    list(U = c(MU, PU), PUW = PUW, MUW = MUW)    
}
## Initializing objects
Draws <- 100
Idx <- 1:Draws
Samples <- c(12, 18, 24, 30, 36, 42, 48, 54, 60)
SS <- length(Samples)
PU <- matrix(NA, ncol = length(Samples), nrow = Draws)
MU <- matrix(NA, ncol = length(Samples), nrow = Draws)
colnames(PU) <- colnames(MU) <- paste("S", Samples, sep = "")
PUW <- array(NA, dim = c(Draws, N, SS))
MUW <- array(NA, dim = c(Draws, N, SS))
set.seed(123456)
RDrawList <- lapply(Idx, function(x) mvrnorm(n = max(Samples),
                                             mu = muSample, Sigma = sigmaSample))
## Carrying out simulation
for(i in 1:length(Samples)){
  cat(paste("Computing for Sample Size", Samples[i], "\n"))
  SampleL <- Samples[i]
  ListData <- lapply(Idx, function(x) RDrawList[[x]][1:Samples[i], ])
  SimResults <- lapply(ListData, PuMuSim)
  MU[, i] <- unlist(lapply(SimResults, function(x) x$U[1]))
  PU[, i] <- unlist(lapply(SimResults, function(x) x$U[2]))
  PUW[, , i] <- matrix(unlist(lapply(SimResults, function(x) x$PUW)),
                       ncol = N, nrow = Draws, byrow = TRUE)
  MUW[, , i] <- matrix(unlist(lapply(SimResults, function(x) x$MUW)),
                       ncol = N, nrow = Draws, byrow = TRUE)
}
