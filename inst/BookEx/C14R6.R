McmcLength <- 5e4
## adaptMCMC
ans1 <- MCMC(p = PUL, n = McmcLength, init = w0, mu = muSample,
             Sigma = sigmaSample, Lambda = riskA, L = targetR,
             iperiod = Nt, nu = convR, acc.rate = 0.5)
ans1$acceptance.rate
w1 <- colMeans(ans1$samples)
w1 <- w1 / sum(w1) * 100
w1
## MCMCpack
ans2 <- MCMCmetrop1R(fun = PUL, theta.init = w0, mcmc = McmcLength, mu = muSample,
                     V = 1e-3 * diag(N), Sigma = sigmaSample, Lambda = riskA, L = targetR,
                     iperiod = Nt, nu = convR)
w2 <- colMeans(ans2)
w2 <- w2 / sum(w2) * 100
w2
## mcmc
ans3 <- metrop(obj = PUL, initial = w0, nbatch = McmcLength, blen = 1,
               scale = 0.025, mu = muSample, Sigma = sigmaSample,
               Lambda = riskA, L = targetR, iperiod = Nt, nu = convR)
ans3$accept
w3 <- colMeans(ans3$batch)
w3 <- w3 / sum(w3) * 100
w3
## rstan
pudat <- list(N = N, Lambda = riskA,  mu = muSample,
              Sigma = sigmaSample, L = targetR,
              iperiod = Nt, nu = convR)
nchain <- 4
StanFile <- file.path(find.package("FRAPO"), "BookEx", "C14S1.stan")
ans4 <- stan(file = StanFile, data = pudat, iter = McmcLength, chains = nchain)
w4 <- drop(get_posterior_mean(ans4, pars = "w")[, nchain + 1]) * 100
w4
## Summarizing results
Wall <- round(rbind(wMEU, w1, w2, w3, w4), 2)
rownames(Wall) <- c("MEU", "adaptMCMC", "MCMCpack", "mcmc", "rstan")
CR <- apply(Wall, 1, function(x) cr(x, sigmaSample))
Wans <- cbind(Wall, CR)
colnames(Wans) <- c("GSPC", "RUA", "GDAXI", "FTSE", "Concentration")
Wans
