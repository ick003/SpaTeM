library(coda)

mcmcChain = sptmRes2mcmc(ResGibbs, c(12,16, 48:50))
gelman.diag(mcmcChain)
gelman.plot(mcmcChain)

mcmcChain = sptmRes2mcmc(ResGibbsM, c(17,33,49,65,93,185:188))
gelman.diag(mcmcChain)
gelman.plot(mcmcChain)

mcmcChain = sptmRes2mcmc(ResGibbsMS, c(187))
gelman.diag(mcmcChain)
gelman.plot(mcmcChain)
