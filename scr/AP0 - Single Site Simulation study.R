library(LearnBayes)
library(mvtnorm)
library(compositions)
library(splines)
library(splines2)
library(slam)
library(Rcpp)
library(deldir)
library(tictoc)
library(tidyverse)
library(mixtools)
library(ggplot2)
library(mgcv)
library(RColorBrewer)


source("./R/functions.R")
files.sources = list.files("R/")
sapply(paste0("R/",files.sources), source)

sourceCpp("./src/cpp_functions.cpp")


set.seed(101)
genData = createSpTdata(max.date = as.Date("31-12-2011", format = "%d-%m-%Y"),
                        min.date = as.Date("01-01-2008", format = "%d-%m-%Y"),
                        parameters = list(sigma = 1*1e-1, tau = 0.1, phi = 1),
                        nSite = 1, nLandUse = 1, byDate = "week", missingMeasures = list(status = T, ratio = 0.25, type= "MAR"))

ggplot(data = genData$df, aes(x = date, y = obs, col = ID)) + geom_point() + 
  geom_segment(aes(x = date, xend = date, y = min(obs), yend = obs), size=0.2) + facet_wrap(~ID)

df = SPTMData(df.obs =genData$df, tempBasis = "bs", tempPeriod = c("%m","%Y"), nSplines = c(12,5),  splinesType = "poly")


coordDF = genData$X[,1:3]
names(coordDF)[1] <- "site"

# Single model 

SpTcov = NULL
df.sptmod = SPTModel(df.sptm = df, coordinates = coordDF, SpTcov = SpTcov)
set.seed(1)

ResGibbs = estimGibbs(df.sptmod, priors = list(beta = list(m0 = 1, s0 = 2, dist = "gauss", constrained = F),
                                               alpha = list(m0 = 0, s0 = 1, dist = "gauss", constrained = F),
                                               gamma = list(m0 = 0, s0 = 2, dist = "gauss"),
                                               sigma2 = list(a0=1, b0=1e-3, dist = "gamma"),
                                               tau = list(a0 = 5, b0= 1, dist = "gamma"),
                                               phi = list(inf = 0.1, sup = 50, dist = "unif"),
                                               pi = list(alpha0= matrix(1, nrow=nlevels(df$obs.data$ID), ncol=1), #
                                                         dist = "dirichlet"),
                                               tauT = list(a0 = 2, b0= 1, dist = "gamma"),
                                               phiT = list(inf = 0.25, sup = 0.5, dist = "unif"),
                                               rho=list(inf=2, sup = 3, dist="unif")),
                      N.run =1500, debug =F,tempRE = "notcorr", print.res = T, 
                      nBatch = 2, parallel = F, nCluster = 1)

keepRun = 800:1500
plotSPTM(ResGibbs, "tempBasis", keepRun = keepRun)
plotSPTM(ResGibbs, "beta")
plotSPTM(ResGibbs, "alpha")
plotSPTM(ResGibbs, "spatialBasis")
plotSPTM(ResGibbs, "hyperparameter", keepRun = keepRun)
plotSPTM(ResGibbs, "residuals", keepRun = keepRun)
plotSPTM(ResGibbs, "covariates", keepRun = keepRun)


predData = genData$df
predData$yP =  predictSPTM(ResGibbs, keepRun = keepRun)$Yp
predData$CI025 = apply(ResGibbs$GibbsOut$yHat[keepRun,,1],2,quantile,0.025)
predData$CI975 = apply(ResGibbs$GibbsOut$yHat[keepRun,,1],2,quantile,0.975)
predData$yPp =  predictSPTM(ResGibbs, keepRun = keepRun)$YpwN

ggplot(data = predData, aes(x = date, y = obs, col = ID)) + geom_point() + geom_line() + 
  geom_line(aes(x = date, y = yP, group = ID), col = "black") +
  geom_line(aes(x = date, y = yPp, group = ID), col = "red") +
  facet_wrap(~ID) + 
  geom_ribbon(aes(ymin = CI025, ymax = CI975), fill = "grey70", col = "grey70", alpha = 0.2)

sum((predData$obs - predData$yP)^2)

colSums(apply(ResGibbs$GibbsOut$yHat[keepRun,,1],1, function(x) (x - ResGibbs$GibbsOut$y$obs)^2))

hist(rigamma(1000, nrow(predData)/2 + 1,sum((predData$obs - predData$yP)^2)/2 + 1e-3 ), freq=F, xlim = c(0,0.02))
hist(ResGibbs$GibbsOut$theta$sigma2H[keepRun,,1], add=T, col = "red", freq = F, breaks = seq(0,5,0.001))

basisSplines = do.call(cbind,ResGibbs$GibbsOut$basis$splines)
basis = ResGibbs$GibbsOut$basis
nSplines = ncol(basisSplines)
nBasis = 2;nComp = 1


dfmp = NULL
for(kC in 1){
  for(i in 1:nBasis){
    t0 = ResGibbs$GibbsOut$y$date
    xPred0 = seq(min(t0), max(t0), "week")
    if(basis$tempPeriod[i] == "%m"){
      xPred0 = as.Date(paste("2014-",format(xPred0,"%m-%d"), sep = ""))
    }
    xBasis = predict(basis$splines[[i]],xPred0)
    yPred = xBasis %*% t(ResGibbs$GibbsOut$theta$alphaH[keepRun,basis$idx[[i]],kC,1])
    xPred = seq(min(ResGibbs$GibbsOut$y$date), max(ResGibbs$GibbsOut$y$date), by = "week")
    dfmp = rbind(dfmp, data.frame(date = xPred, f = rowMeans(yPred), 
                                  f025 = apply(yPred,1,quantile,0.05), f975 = apply(yPred,1,quantile,0.95), 
                                  basis = rep(basis$tempPeriod[i], length(xPred)),
                                  landuse = rep("Cluster 1", length(xPred))))
  }
}



ggplot(data = dfmp) + geom_line(aes(x = date, y = f, group = basis, color = basis)) + 
  geom_ribbon(aes(x = date, y = f, ymin = f025, ymax = f975, group = basis, fill = basis), alpha = 0.3)+
  facet_grid(basis~landuse) + theme_bw() + scale_color_grey() + scale_fill_grey() + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = 'none')


