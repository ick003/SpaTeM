library(LearnBayes)
library(mvtnorm)
library(compositions)
library(splines)
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



genData = createSpTdata(max.date = as.Date("31-12-2011", format = "%d-%m-%Y"),
                        min.date = as.Date("01-01-2008", format = "%d-%m-%Y"),
                        parameters = list(sigma = 1*1e-1, tau = 0.1, phi = 1),
                        nSite = 9, nLandUse = 3, byDate = "month", missingMeasures = list(status = T, ratio = 0.1, type= "MAR"))

ggplot(data = genData$df, aes(x = date, y = obs, col = ID)) + geom_point() + 
  geom_segment(aes(x = date, xend = date, y = 0, yend = obs), size=0.2) + facet_wrap(~ID)

df = SPTMData(df.obs =genData$df, tempBasis = "bs", tempPeriod = c("%m","%Y"), nSplines = c(11,5),  splinesType = "poly")

coordDF = genData$X[,1:3]
names(coordDF)[1] <- "site"

# Single model 

CLR = clr(genData$luComp)

SpTcov = cbind(CLR[match(df$obs.data$ID,rownames(genData$luComp)),],1)
df.sptmod = SPTModel(df.sptm = df, coordinates = coordDF, SpTcov = SpTcov)
set.seed(1)

ResGibbs = estimGibbs(df.sptmod, priors = list(beta = list(m0 = 1, s0 = 1, dist = "gauss", constrained = F),
                                               alpha = list(m0 = 0, s0 = 1, dist = "gauss", constrained = T),
                                               gamma = list(m0 = 0, s0 = 1, dist = "gauss"),
                                               sigma2 = list(a0=20, b0=8, dist = "gamma"),
                                               tau = list(a0 = 5, b0= 1, dist = "gamma"),
                                               phi = list(inf = 0.1, sup = 50, dist = "unif"),
                                               pi = list(alpha0= matrix(1, nrow=nlevels(df$obs.data$ID), ncol=1), #
                                                         dist = "dirichlet"),
                                               tauT = list(a0 = 2, b0= 1, dist = "gamma"),
                                               phiT = list(inf = 0.25, sup = 0.5, dist = "unif"),
                                               rho=list(inf=2, sup = 3, dist="unif")),
                      N.run = 500, debug =F,tempRE = "notcorr", print.res = T, 
                      nBatch = 2, parallel = F, nCluster = 1)

plotSPTM(ResGibbs, "tempBasis")
plotSPTM(ResGibbs, "beta")
plotSPTM(ResGibbs, "alpha")
plotSPTM(ResGibbs, "spatialBasis")
plotSPTM(ResGibbs, "hyperparameter")
plotSPTM(ResGibbs, "residuals")
plotSPTM(ResGibbs, "covariates")

keepRun = 200:500
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

# Mixture model

df.sptmod = SPTModel(df.sptm = df, coordinates = coordDF, SpTcov = NULL)
pi0 = genData$luComp[,1:3]
pi0 = pi0/ rowSums(pi0)#genData$grTrue#
set.seed(1)
ResGibbsM = estimGibbs(df.sptmod, priors = list(beta = list(m0 = 1, s0 = 1, dist = "gauss", constrained = F),
                                                alpha = list(m0 = 0, s0 = 1, dist = "gauss", constrained = T),
                                                gamma = list(m0 = 0, s0 = 1, dist = "gauss"),
                                                sigma2 = list(a0=20, b0=8, dist = "igamma"),
                                                tau = list(a0 = 5, b0= 1, dist = "igamma"),
                                                phi = list(inf = 0.1, sup = 50, dist = "unif"),
                                                tauT = list(a0 = 2, b0= 1, dist = "gamma"),
                                                phiT = list(inf = 0.1, sup = 0.8, dist = "unif"),
                                                rho = list(inf = 0.1, sup = 0.8, dist = "unif"),
                                                pi = list(alpha0 = pi0, dist = "dirichlet")),
                       N.run = 500, debug = F, tempRE = "notcorr", model = "simpleMixture",
                       print.res = T, nBatch = 2,
                       parallel = F, nCluster = ncol(pi0))


plotSPTM(ResGibbsM, "tempBasis", keepRun = 200:500)
plotSPTM(ResGibbsM, "beta")
plotSPTM(ResGibbsM, "alpha")
plotSPTM(ResGibbsM, "spatialBasis")
plotSPTM(ResGibbsM, "hyperparameter", keepRun = 250:500)
plotSPTM(ResGibbsM, "residuals", keepRun = 250:500)
plotSPTM(ResGibbsM, "covariates", keepRun = 250:500)
plotSPTM(ResGibbsM, "cluster", keepRun = 1250:2500)
plotSPTM(ResGibbsM, "map")



apply(t(apply(ResGibbsM$GibbsOut$theta$zH[,,,1],c(2,3),mean)),1,which.max)
genData$grTrue
genData$betaTrue

ResGibbsM$GibbsOut$logPostDist$ll[250,,]

predData = genData$df
predData$yP =  predictSPTM(ResGibbsM, keepRun = 250:500)$Yp
predData$CI025 = apply(ResGibbsM$GibbsOut$yHat[250:500,,1],2,quantile,0.025)
predData$CI975 = apply(ResGibbsM$GibbsOut$yHat[250:500,,1],2,quantile,0.975)
predData$yPp =  predictSPTM(ResGibbsM, keepRun = 250:500)$YpwN

Cldf = data.frame(ID = 1:9,Cluster = apply(t(apply(ResGibbsM$GibbsOut$theta$piH[250:500,,,2],c(2,3),mean)),1,which.max))

predData = merge(predData, Cldf, by  ="ID", all.x = TRUE)


ggplot(data = predData[predData$Cluster %in% c(1,2,3),], aes(x = date, y = obs, col = ID)) + geom_point() + geom_line() + 
  geom_line(aes(x = date, y = yP, group = ID), col = "black") +
  geom_line(aes(x = date, y = yPp, group = ID), col = "red") +
  facet_wrap(Cluster~ID) + 
  geom_ribbon(aes(ymin = CI025, ymax = CI975), fill = "grey70", col = "grey70", alpha = 0.2)


mean(exp(ResGibbs$GibbsOut$logPostDist$ll[keepRun,,2]))/mean(exp(ResGibbsM$GibbsOut$logPostDist$ll[keepRun,,2]))

# Spatial mixture 


dd = deldir(genData$coordinates[,1], genData$coordinates[,2],suppressMsge=TRUE)
tile.dd = tile.list(dd)
Bds = comPts(tile.dd)


Blocks = list(c(1,2,6),c(3,5,7,8), c(4,9))
set.seed(1)
ResGibbsMS = estimGibbs(df.sptmod, priors = list(beta = list(m0 = 1, s0 = 0.1, dist = "gauss", constrained =F),
                                                 alpha = list(m0 = 0, s0 = 2, dist = "gauss", constrained = T),# converged with F
                                                 gamma = list(m0 = 0, s0 = 0.1, dist = "gauss"),
                                                 sigma2 = list(a0=2, b0=1, dist = "gamma"),
                                                 tau = list(a0 = 5, b0= 1, dist = "gamma"),
                                                 phi = list(inf = 0.1, sup = 50, dist = "unif"),
                                                 tauT = list(a0 = 2, b0= 1, dist = "gamma"),
                                                 phiT = list(inf = 0.1, sup = 0.8, dist = "unif"),
                                                 pi = list(alpha0 = pi0, dist = "dirichlet"),
                                                 rho=list(inf=0, sup = 3, dist="unif")),
                        N.run = 500, debug = F, model = "spatialMixture",tempRE = "notcorr", 
                        print.res = T, nBatch = 2, Blocks = Blocks,
                        parallel = F, nCluster = ncol(pi0))

plotSPTM(ResGibbsMS, "tempBasis", keepRun = 200:500)
plotSPTM(ResGibbsMS, "beta")
plotSPTM(ResGibbsMS, "alpha")
plotSPTM(ResGibbsMS, "spatialBasis")
plotSPTM(ResGibbsMS, "hyperparameter", keepRun = 250:500)
plotSPTM(ResGibbsMS, "residuals", keepRun = 250:500)
plotSPTM(ResGibbsMS, "covariates", keepRun = 250:500)
plotSPTM(ResGibbsMS, "cluster", keepRun = 250:500)
plotSPTM(ResGibbsMS, "map")



apply(t(apply(ResGibbsMS$GibbsOut$theta$zH[,,,1],c(2,3),mean)),1,which.max)
genData$grTrue
genData$betaTrue

ResGibbsM$GibbsOut$logPostDist$ll[500,,]

predData = genData$df
predData$yP =  predictSPTM(ResGibbsMS, keepRun = 250:500)$Yp
predData$CI025 = apply(ResGibbsMS$GibbsOut$yHat[250:500,,1],2,quantile,0.025)
predData$CI975 = apply(ResGibbsMS$GibbsOut$yHat[250:500,,1],2,quantile,0.975)
predData$yPp =  predictSPTM(ResGibbsMS, keepRun = 250:500)$YpwN

Cldf = data.frame(ID = 1:9,Cluster = apply(t(apply(ResGibbsMS$GibbsOut$theta$piH[250:500,,,2],c(2,3),mean)),1,which.max))

predData = merge(predData, Cldf, by  ="ID", all.x = TRUE)


ggplot(data = predData[predData$Cluster %in% c(1,2,3),], aes(x = date, y = obs, col = ID)) + geom_point() + geom_line() + 
  geom_line(aes(x = date, y = yP, group = ID), col = "black") +
  geom_line(aes(x = date, y = yPp, group = ID), col = "red") +
  facet_wrap(Cluster~ID) + 
  geom_ribbon(aes(ymin = CI025, ymax = CI975), fill = "grey70", col = "grey70", alpha = 0.2)


mean(exp(ResGibbsMS$GibbsOut$logPostDist$ll[keepRun,,2]))/mean(exp(ResGibbsM$GibbsOut$logPostDist$ll[keepRun,,2]))
