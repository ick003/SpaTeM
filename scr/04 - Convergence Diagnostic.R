library(coda)

mcmcChain = sptmRes2mcmc(ResGibbs)
gelman.diag(mcmcChain)
gelman.plot(mcmcChain)

mcmcChain = sptmRes2mcmc(ResGibbsM)
gelman.diag(mcmcChain)
gelman.plot(mcmcChain)

mcmcChain = sptmRes2mcmc(ResGibbsMS)
gelman.diag(mcmcChain)
gelman.plot(mcmcChain)