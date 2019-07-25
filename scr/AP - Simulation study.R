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


set.seed(102)
genData = createSpTdata(max.date = as.Date("31-12-2011", format = "%d-%m-%Y"),
                        min.date = as.Date("01-01-2008", format = "%d-%m-%Y"),
                        parameters = list(sigma = 1*1e-2, tau = 0.1, phi = 1),
                        nSite = 9, nLandUse = 3, byDate = "month", missingMeasures = list(status = T, ratio = 0.25, type= "MAR"))
genData$grTrue

ggplot(data = genData$df, aes(x = date, y = obs, col = ID)) + geom_point() + 
  geom_segment(aes(x = date, xend = date, y = min(obs), yend = obs), size=0.2) + facet_wrap(~ID)

df = SPTMData(df.obs =genData$df, tempBasis = "bs", tempPeriod = c("%m","%Y"), nSplines = c(12,5),  splinesType = "poly")

coordDF = genData$X[,1:3]
names(coordDF)[1] <- "site"

# Single model 

CLR = clr(genData$luComp[,-4])

SpTcov = cbind(CLR[match(df$obs.data$ID,rownames(genData$luComp)),],1)
df.sptmod = SPTModel(df.sptm = df, coordinates = coordDF, SpTcov = SpTcov)
set.seed(1)

ResGibbs = estimGibbs(df.sptmod, priors = list(beta = list(m0 = 1, s0 = 2, dist = "gauss", constrained = F),
                                               alpha = list(m0 = 0, s0 = 1, dist = "gauss", constrained = F),
                                               gamma = list(m0 = 0, s0 = 2, dist = "gauss"),
                                               sigma2 = list(a0=20, b0=8, dist = "gamma"),
                                               tau = list(a0 = 5, b0= 1, dist = "gamma"),
                                               phi = list(inf = 0.1, sup = 50, dist = "unif"),
                                               pi = list(alpha0= matrix(1, nrow=nlevels(df$obs.data$ID), ncol=1), #
                                                         dist = "dirichlet"),
                                               tauT = list(a0 = 2, b0= 1, dist = "gamma"),
                                               phiT = list(inf = 0.25, sup = 0.5, dist = "unif"),
                                               rho=list(inf=2, sup = 3, dist="unif")),
                      N.run = 1500, debug =F,tempRE = "notcorr", print.res = T, 
                      nBatch = 2, parallel = F, nCluster = 1)

Nrun = 1500
idxBatch = which.max(ResGibbs$GibbsOut$logPostDist$ll[Nrun,,])
keepRun = round(0.5*Nrun):Nrun
plotSPTM(ResGibbs, "tempBasis", keepRun = keepRun)
plotSPTM(ResGibbs, "beta")
plotSPTM(ResGibbs, "alpha")
plotSPTM(ResGibbs, "spatialBasis")
plotSPTM(ResGibbs, "hyperparameter", keepRun = keepRun)
plotSPTM(ResGibbs, "residuals", keepRun = keepRun)
plotSPTM(ResGibbs, "covariates", keepRun = keepRun)


colMeans(ResGibbs$GibbsOut$theta$betaIH[keepRun,,1])

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
ResGibbsM = estimGibbs(df.sptmod, priors = list(beta = list(m0 = 1, s0 = 2, dist = "gauss", constrained = F),
                                                alpha = list(m0 = 0, s0 = 2, dist = "gauss", constrained = F),
                                                gamma = list(m0 = 0, s0 = 2, dist = "gauss"),
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

Nrun = 500
idxBatchM = which.max(ResGibbsM$GibbsOut$logPostDist$ll[Nrun,,])
Prdf = apply(ResGibbsM$GibbsOut$theta$zProbH[round(0.5*Nrun):Nrun,,,idxBatchM],c(2,3),mean)
Cldf = data.frame(ID = 1:9,Cluster = apply(t(apply(ResGibbsM$GibbsOut$theta$zH[round(0.5*Nrun):Nrun,,,idxBatchM],c(2,3),mean)),1,which.max))
idxKeepM = which(colSums(apply(ResGibbsM$GibbsOut$theta$zH[,,,idxBatchM], 1, function(x) apply(x, 2,which.max) - Cldf$Cluster))==0)
idxKeepM = idxKeepM[idxKeepM > round(0.5*Nrun)]

plotSPTM(ResGibbsM, "tempBasis", keepRun = idxKeepM)
plotSPTM(ResGibbsM, "beta")
plotSPTM(ResGibbsM, "alpha")
plotSPTM(ResGibbsM, "spatialBasis")
plotSPTM(ResGibbsM, "hyperparameter", keepRun = idxKeepM)
plotSPTM(ResGibbsM, "residuals", keepRun = idxKeepM)
plotSPTM(ResGibbsM, "covariates", keepRun = idxKeepM)
plotSPTM(ResGibbsM, "cluster", keepRun = idxKeepM)

predData = genData$df
predData$yP =  predictSPTM(ResGibbsM, keepRun = idxKeepM)$Yp
predData$CI025 = apply(ResGibbsM$GibbsOut$yHat[idxKeepM,,idxBatchM],2,quantile,0.025)
predData$CI975 = apply(ResGibbsM$GibbsOut$yHat[idxKeepM,,idxBatchM],2,quantile,0.975)
predData$yPp =   apply(ResGibbsM$GibbsOut$yHat[idxKeepM,,idxBatchM],2,mean)

predData = merge(predData, Cldf, by  ="ID", all.x = TRUE)

ggplot(data = predData[predData$Cluster %in% c(1,2,3),], aes(x = date, y = obs, col = ID)) + geom_point() + geom_line() + 
  geom_line(aes(x = date, y = yP, group = ID), col = "black") +
  geom_line(aes(x = date, y = yPp, group = ID), col = "red") +
  facet_wrap(Cluster~ID) + 
  geom_ribbon(aes(ymin = CI025, ymax = CI975), fill = "grey70", col = "grey70", alpha = 0.2)


# Spatial mixture 


dd = deldir(genData$coordinates[,1], genData$coordinates[,2],suppressMsge=TRUE)
tile.dd = tile.list(dd)
Bds = comPts(tile.dd)

par(mfrow=c(1,1))
plot(genData$coordinates, cex=0.1)
text(genData$coordinates, labels = 1:9)
Bds
Blocks = list(c(1,2,7),c(3,5,9), c(6,8),c(4))

set.seed(1)
ResGibbsMS = estimGibbs(df.sptmod, priors = list(beta = list(m0 = 1, s0 = 2, dist = "gauss", constrained =T),
                                                 alpha = list(m0 = 0, s0 = 2, dist = "gauss", constrained = T),# converged with F
                                                 gamma = list(m0 = 0, s0 = 2, dist = "gauss"),
                                                 sigma2 = list(a0=2, b0=1, dist = "gamma"),
                                                 tau = list(a0 = 5, b0= 1, dist = "gamma"),
                                                 phi = list(inf = 0.1, sup = 50, dist = "unif"),
                                                 tauT = list(a0 = 2, b0= 1, dist = "gamma"),
                                                 phiT = list(inf = 0.1, sup = 0.8, dist = "unif"),
                                                 pi = list(alpha0 = pi0, dist = "dirichlet"),
                                                 rho=list(inf=0, sup = 3, dist="unif")),
                        N.run = 2500, debug = F, model = "spatialMixture",tempRE = "notcorr", 
                        print.res = T, nBatch = 2, Blocks = Blocks,
                        parallel = F, nCluster = ncol(pi0))

Nrun = 2500
idxBatchMS = which.max(ResGibbsMS$GibbsOut$logPostDist$ll[Nrun,,])
Cldf = data.frame(ID = 1:9,Cluster = apply(t(apply(ResGibbsMS$GibbsOut$theta$zH[round(0.5*Nrun):Nrun,,,idxBatchMS],c(2,3),mean)),1,which.max))
idxKeepMS = which(colSums(apply(ResGibbsMS$GibbsOut$theta$zH[,,,idxBatchMS], 1, function(x) apply(x, 2,which.max) - Cldf$Cluster))==0)
idxKeepMS = idxKeepMS[idxKeepMS > round(0.5*Nrun)]

plotSPTM(ResGibbsMS, "tempBasis", keepRun = idxKeepMS)
plotSPTM(ResGibbsMS, "beta")
plotSPTM(ResGibbsMS, "alpha")
plotSPTM(ResGibbsMS, "spatialBasis")
plotSPTM(ResGibbsMS, "hyperparameter", keepRun = idxKeepMS)
plotSPTM(ResGibbsMS, "residuals", keepRun = idxKeepMS)
plotSPTM(ResGibbsMS, "covariates", keepRun =idxKeepMS)
plotSPTM(ResGibbsMS, "cluster", keepRun = idxKeepMS)

predData = genData$df
predData$yP =  predictSPTM(ResGibbsMS, keepRun = idxKeepMS)$Yp
predData$CI025 = apply(ResGibbsMS$GibbsOut$yHat[idxKeepMS,,idxBatchMS],2,quantile,0.025)
predData$CI975 = apply(ResGibbsMS$GibbsOut$yHat[idxKeepMS,,idxBatchMS],2,quantile,0.975)
predData$yPp = apply(ResGibbsMS$GibbsOut$yHat[idxKeepMS,,idxBatchMS],2,mean)

predData = merge(predData, Cldf, by  ="ID", all.x = TRUE)


ggplot(data = predData[predData$Cluster %in% c(1,2,3),], aes(x = date, y = obs, col = ID)) + geom_point() + geom_line() + 
  geom_line(aes(x = date, y = yP, group = ID), col = "black") +
  geom_line(aes(x = date, y = yPp, group = ID), col = "red") +
  facet_wrap(Cluster~ID) + 
  geom_ribbon(aes(ymin = CI025, ymax = CI975), fill = "grey70", col = "grey70", alpha = 0.2)


#### MAP Estimates ###

library(reshape2)

nRuns = dim(ResGibbs$GibbsOut$theta$alphaH)[1]
nChains = dim(ResGibbs$GibbsOut$theta$alphaH)[4]

PointEst = as.data.frame(t(round(rbind(apply(ResGibbs$GibbsOut$theta$gammaH[keepRun,,idxBatch], 2, mean),
                                       apply(ResGibbs$GibbsOut$theta$gammaH[keepRun,,idxBatch], 2, quantile, c(0.025, 0.975))),3)))

colnames(PointEst)[1] <- "Post. mean"
PointEst$LandUse <- c(colnames(CLR), "Intercept")

PointEst$Model <- "CLR"

ggplot(data = PointEst) + geom_pointrange(aes(x = reorder(LandUse, `Post. mean`), y = `Post. mean`, ymin = `2.5%`, ymax = `97.5%`))


nRuns = dim(ResGibbsM$GibbsOut$theta$alphaH)[1]
nChains = dim(ResGibbsM$GibbsOut$theta$alphaH)[4]

PointEstM = as.data.frame(t(round(rbind(apply(ResGibbsM$GibbsOut$theta$gammaH[idxKeepM,,idxBatchM], 2, mean),
                                        apply(ResGibbsM$GibbsOut$theta$gammaH[idxKeepM,,idxBatchM], 2, quantile, c(0.025, 0.975))),3)))

colnames(PointEstM)[1] <- "Post. mean"
PointEstM$LandUse <- c(colnames(CLR))

PointEstM$Model <- "Mixture"

nRuns = dim(ResGibbsMS$GibbsOut$theta$alphaH)[1]
nChains = dim(ResGibbsMS$GibbsOut$theta$alphaH)[4]

PointEstMS = as.data.frame(t(round(rbind(apply(ResGibbsMS$GibbsOut$theta$gammaH[idxKeepMS,,idxBatchMS], 2, mean),
                                         apply(ResGibbsMS$GibbsOut$theta$gammaH[idxKeepMS,,idxBatchMS], 2, quantile, c(0.025, 0.975))),3)))

colnames(PointEstMS)[1] <- "Post. mean"
PointEstMS$LandUse <- c(colnames(CLR))

PointEstMS$Model <- "Spatial Mixture"

PE = rbind(PointEst[-4,], PointEstM, PointEstMS)


ggplot(data = PE) + geom_pointrange(aes(x = reorder(LandUse, `Post. mean`), 
                                        y = `Post. mean`, ymin = `2.5%`, ymax = `97.5%`, shape = Model, linetype = Model), 
                                    position=position_dodge(width=c(0.4)), size= 0.5) + theme_bw()+
  theme(axis.text.x = element_text(angle = 40, hjust = 0.6, vjust = 0.6, size=10),axis.title.x=element_blank(), axis.title.y=element_blank())


# Temporal functions

basisSplines = do.call(cbind,ResGibbsMS$GibbsOut$basis$splines)
basis = ResGibbsMS$GibbsOut$basis
nSplines = ncol(basisSplines)
nBasis = 2;nComp = 3


dfmp = NULL
for(kC in c(1,2,3)){
  for(i in 1:nBasis){
    t0 = ResGibbsMS$GibbsOut$y$date
    xPred0 = seq(min(t0), max(t0), "month")
    if(basis$tempPeriod[i] == "%m"){
      xPred0 = as.Date(paste("2014-",format(xPred0,"%m-%d"), sep = ""))
    }
    xBasis = predict(basis$splines[[i]],xPred0)
    yPred = xBasis %*% t(ResGibbsMS$GibbsOut$theta$alphaH[idxKeepMS,basis$idx[[i]],kC,idxBatchMS])
    xPred = seq(min(ResGibbsMS$GibbsOut$y$date), max(ResGibbsMS$GibbsOut$y$date), by = "month")
    dfmp = rbind(dfmp, data.frame(date = xPred, f = rowMeans(yPred), 
                                  f025 = apply(yPred,1,quantile,0.025), f975 = apply(yPred,1,quantile,0.975), 
                                  basis = rep(basis$tempPeriod[i], length(xPred)),
                                  landuse = rep(names(CLR)[kC], length(xPred))))
  }
}



ggplot(data = dfmp) + geom_line(aes(x = date, y = f, group = basis)) + 
  geom_ribbon(aes(x = date, y = f, ymin = f025, ymax = f975, group = basis), alpha = 0.3)+
  facet_grid(basis~landuse) + theme_bw() + scale_color_grey() + scale_fill_grey() + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = 'none')

