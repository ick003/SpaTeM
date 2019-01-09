library(coda)

mcmcChain = sptmRes2mcmc(ResGibbs,67:71)
gelman.diag(mcmcChain)
gelman.plot(mcmcChain)

mcmcChain = sptmRes2mcmc(ResGibbsM, 329:330)
gelman.diag(mcmcChain)
gelman.plot(mcmcChain)

mcmcChain = sptmRes2mcmc(ResGibbsMS, c(187))
gelman.diag(mcmcChain)
gelman.plot(mcmcChain)
