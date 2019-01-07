library(coda)

mcmcChain = sptmRes2mcmc(ResGibbs, c(62:68))
gelman.diag(mcmcChain)
gelman.plot(mcmcChain)

mcmcChain = sptmRes2mcmc(ResGibbsM)
gelman.diag(mcmcChain)
gelman.plot(mcmcChain)

mcmcChain = sptmRes2mcmc(ResGibbsMS, c(187))
gelman.diag(mcmcChain)
gelman.plot(mcmcChain)
