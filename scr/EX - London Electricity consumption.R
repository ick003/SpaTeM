library(OpenML)
library(farff)
library(TSrepr)
library(reshape2)


data <- OpenML::getOMLDataSet(data.id = 41060)
data <- data.matrix(data$data)
nID = 8
data_cons <- data[1:nID,]

data_ave_prof <- repr_matrix(data_cons,
                             func = repr_seas_profile,
                             args = list(freq = 48,
                                         func = median),
                             normalise = TRUE,
                             func_norm = norm_z)

colnames(data_ave_prof) <- seq(0,23.5,0.5)
rownames(data_ave_prof) <- paste("Consumer", 1:nID)
data_plot <- (melt(data_ave_prof))
colnames(data_plot) <- c("ID", "date", "obs")

ggplot(data_plot, aes(x = date, y = obs)) + geom_point(aes(group = ID), col = "grey30", alpha = 0.3)+
  stat_smooth(span = 0.85, col = "red") + facet_wrap(~ID)




df = SPTMData(df.obs =data_plot, tempBasis = "bs", tempPeriod = c("%H"), nSplines = c(24),  splinesType = "poly")

coordDF = data.frame(site = unique(data_plot$ID),
                     longitude = runif(nID,0,1),
                     latitude = runif(nID,0,1))

#### Spline based model ####

# Single model 

SpTcov = NULL
df.sptmod = SPTModel(df.sptm = df, coordinates = coordDF, SpTcov = SpTcov)
set.seed(1)

ResGibbs = estimGibbs(df.sptmod, priors = list(beta = list(m0 = 1, s0 = 0.2, dist = "gauss", constrained = F),
                                               alpha = list(m0 = 0, s0 = 0.1, dist = "gauss", constrained = F),
                                               gamma = list(m0 = 0, s0 = 0.2, dist = "gauss"),
                                               sigma2 = list(a0=0.8, b0=2, dist = "gamma"),
                                               tau = list(a0 = 5, b0= 1, dist = "gamma"),
                                               phi = list(inf = 0.1, sup = 50, dist = "unif"),
                                               pi = list(alpha0= matrix(1, nrow=nlevels(df$obs.data$ID), ncol=1), #
                                                         dist = "dirichlet"),
                                               tauT = list(a0 = 2, b0= 1, dist = "gamma"),
                                               phiT = list(inf = 0.25, sup = 0.5, dist = "unif"),
                                               rho=list(inf=2, sup = 3, dist="unif")),
                      N.run =500, debug =F,tempRE = "notcorr", print.res = T, 
                      nBatch = 2, parallel = F, nCluster = 1)

keepRun = 300:500
plotSPTM(ResGibbs, "tempBasis", keepRun = keepRun)
plotSPTM(ResGibbs, "beta")
plotSPTM(ResGibbs, "alpha")
plotSPTM(ResGibbs, "spatialBasis")
plotSPTM(ResGibbs, "hyperparameter", keepRun = keepRun)
plotSPTM(ResGibbs, "residuals", keepRun = keepRun)
plotSPTM(ResGibbs, "covariates", keepRun = keepRun)


predData = data_plot
predData$yP =  predictSPTM(ResGibbs, keepRun = keepRun)$Yp
predData$CI025 = apply(ResGibbs$GibbsOut$yHat[keepRun,,1],2,quantile,0.025)
predData$CI975 = apply(ResGibbs$GibbsOut$yHat[keepRun,,1],2,quantile,0.975)
predData$yPp =  predictSPTM(ResGibbs, keepRun = keepRun)$YpwN

ggplot(data = predData, aes(x = date, y = obs, col = ID)) + geom_point() + geom_line() + 
  geom_line(aes(x = date, y = yP, group = ID), col = "black") +
  geom_line(aes(x = date, y = yPp, group = ID), col = "red") +
  facet_wrap(~ID) + 
  geom_ribbon(aes(ymin = CI025, ymax = CI975), fill = "grey70", col = "grey70", alpha = 0.2)

# Multiple CLusters
nCl = 8
pi0 = matrix(runif(8*nCl,0,1), nrow = 8, ncol = nCl)
pi0 = pi0/ rowSums(pi0)#genData$grTrue#
colnames(pi0) <- paste0("Cl", 1:nCl)
rownames(pi0) <- coordDF$site
set.seed(1)
ResGibbsM = estimGibbs(df.sptmod, priors = list(beta = list(m0 = 1, s0 = 2, dist = "gauss", constrained = F),
                                                alpha = list(m0 = 0, s0 = 2, dist = "gauss", constrained = F),
                                                gamma = list(m0 = 0, s0 = 2, dist = "gauss"),
                                                sigma2 = list(a0=200, b0=8, dist = "igamma"),
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
Cldf = data.frame(ID = coordDF$site,Cluster = apply(t(apply(ResGibbsM$GibbsOut$theta$zH[round(0.5*Nrun):Nrun,,,idxBatchM],c(2,3),mean)),1,which.max))
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


predData = data_plot
predData$yP =  predictSPTM(ResGibbsM, keepRun = idxKeepM)$Yp
predData$CI025 = apply(ResGibbsM$GibbsOut$yHat[idxKeepM,,idxBatchM],2,quantile,0.025)
predData$CI975 = apply(ResGibbsM$GibbsOut$yHat[idxKeepM,,idxBatchM],2,quantile,0.975)
predData$yPp =   apply(ResGibbsM$GibbsOut$yHat[idxKeepM,,idxBatchM],2,quantile, 0.5)
predData$yPp =   apply(ResGibbsM$GibbsOut$yHat[idxKeepM,,idxBatchM],2,mean)

predData = merge(predData, Cldf, by  ="ID", all.x = TRUE)

ggplot(data = predData[predData$Cluster %in% unique(Cldf$Cluster)[1:6],], aes(x = date, y = obs, col = ID)) + geom_point() + geom_line() + 
  geom_line(aes(x = date, y = yP, group = ID), col = "black") +
  geom_line(aes(x = date, y = yPp, group = ID), col = "red") +
  facet_wrap(Cluster~ID) + 
  geom_ribbon(aes(ymin = CI025, ymax = CI975), fill = "grey70", col = "grey70", alpha = 0.2)
  

plot(do.call(cbind,ResGibbsM$GibbsOut$basis$splines)[,ResGibbsM$GibbsOut$basis$idx[[1]]]%*% 
  ResGibbsM$GibbsOut$theta$alphaH[Nrun-10,ResGibbsM$GibbsOut$basis$idx[[1]], 1,1], type="l")

#### GP based model ####

df.sptmod = SPTModel(df.sptm = df, coordinates = coordDF, SpTcov = NULL)

# SImple model 

thetaGP_HP = list(k1=list(a1 = 10, b1=20, a2=10, b2 = 4)) 
#                  k2 = list(a3=10, b3=10, a4=1000, b4 = 1,as=10, bs=10, af=10, bf = 3600 )) # q3 = 1, q4 = 1000, qs = 1000, f = 1/360
kernelList = list(k_longterm)
attr(kernelList, "name") <- c("long term")
attr(kernelList, "parameters") <- list(c("q1", "q2"))
attr(kernelList, "type") <- c("temporal")

ResGibbs = estimGibbs(df.sptmod, priors = list(beta = list(m0 = 1, s0 = 0.1, dist = "gauss", constrained =T),
                                               thetaGP = thetaGP_HP,
                                               gamma = list(m0 = 0, s0 = 0.2, dist = "gauss"),
                                               sigma2 = list(a0=100, b0=1, dist = "igamma"),
                                               tau = list(a0 = 100, b0= 1, dist = "igamma"),
                                               phi = list(inf = 0.8, sup = 1.2, dist = "unif"),
                                               pi = list(alpha0= matrix(1, nrow=nlevels(df$obs.data$ID), ncol=1), #
                                                         dist = "dirichlet"),
                                               tauT = list(a0 = 2, b0= 1, dist = "gamma"),
                                               phiT = list(inf = 0.1, sup = 0.8, dist = "unif"),
                                               rho = list(inf = 0.1, sup = 0.8, dist = "unif")),
                      N.run = 500, debug = F, tempRE = "gp", model = "simpleMixture",
                      print.res = T, nBatch = 2,
                      parallel = F, nCluster = 1,kernelList = kernelList)


keepRun = 250:500
plotSPTM(ResGibbs, "spatialBasis")
plotSPTM(ResGibbs, "hyperparameter", keepRun = keepRun)
plotSPTM(ResGibbs, "residuals", keepRun = keepRun)
plotSPTM(ResGibbs, "covariates", keepRun = keepRun)


predData = data_plot
predData$yP =  predictSPTM(ResGibbs, keepRun = keepRun)$Yp
predData$CI025 = apply(ResGibbs$GibbsOut$yHat[keepRun,,1],2,quantile,0.025)
predData$CI975 = apply(ResGibbs$GibbsOut$yHat[keepRun,,1],2,quantile,0.975)
predData$yPp =  predictSPTM(ResGibbs, keepRun = keepRun)$YpwN
predData$CI025p = apply(predictSPTM(ResGibbs, keepRun = keepRun, posterior = T)$YpwN,2,quantile,0.025)
predData$CI975p = apply(predictSPTM(ResGibbs, keepRun = keepRun, posterior = T)$YpwN,2,quantile,0.975)

ggplot(data = predData, aes(x = date, y = obs, col = ID)) + geom_point() + geom_line() + 
  geom_line(aes(x = date, y = yP, group = ID), col = "black") +
  geom_line(aes(x = date, y = yPp, group = ID), col = "red") +
  facet_wrap(~ID) + 
  geom_ribbon(aes(ymin = CI025, ymax = CI975), fill = "grey70", col = "grey50", alpha = 0.2) +
  geom_ribbon(aes(ymin = CI025p, ymax = CI975p), fill = "grey70", col = "grey70", alpha = 0.2)



# Multiple clusters 
nCl = 6
pi0 = matrix(runif(nID*nCl,0,1), nrow = nID, ncol = nCl)
pi0 = pi0/ rowSums(pi0)#genData$grTrue#
colnames(pi0) <- paste0("Cl", 1:nCl)
rownames(pi0) <- coordDF$site
set.seed(1)

thetaGP_HP = list(k1=list(a1 = 10, b1=20, a2=60, b2 = 10)) # q3 = 1, q4 = 1e6, qs = 1e-1, f = 1/360
#plot(0:1440, k_seasonal(0, 0:1440, param = list(q3=1, q4 = 1e6, qs = 1e-1, f = 1/360)))
kernelList = list(k_longterm)
attr(kernelList, "name") <- c("long term")
attr(kernelList, "parameters") <- list(c("q1", "q2"))
attr(kernelList, "type") <- c("temporal")

ResGibbsM = estimGibbs(df.sptmod, priors = list(beta = list(m0 = 1, s0 = 0.1, dist = "gauss", constrained =T),
                                               thetaGP = thetaGP_HP,
                                               gamma = list(m0 = 0, s0 = 0.2, dist = "gauss"),
                                               sigma2 = list(a0=20, b0=1, dist = "igamma"),
                                               tau = list(a0 = 100, b0= 1, dist = "igamma"),
                                               phi = list(inf = 0.8, sup = 1.2, dist = "unif"),
                                               pi = list(alpha0= pi0, dist = "dirichlet"),
                                               tauT = list(a0 = 2, b0= 1, dist = "gamma"),
                                               phiT = list(inf = 0.1, sup = 0.8, dist = "unif"),
                                               rho = list(inf = 0.1, sup = 0.8, dist = "unif")),
                      N.run = 500, debug = F, tempRE = "gp", model = "simpleMixture",
                      print.res = F, nBatch = 2,
                      parallel = F, nCluster = ncol(pi0),kernelList = kernelList)


Nrun = 500
idxBatch = which.max(ResGibbsM$GibbsOut$logPostDist$ll[Nrun,,])
Cldf = data.frame(ID = coordDF$site,Cluster = apply(t(apply(ResGibbsM$GibbsOut$theta$zH[round(0.5*Nrun):Nrun,,,idxBatch],c(2,3),mean)),1,which.max))
idxKeep = which(colSums(apply(ResGibbsM$GibbsOut$theta$zH[,,,idxBatch], 1, function(x) apply(x, 2,which.max) - Cldf$Cluster))==0)
idxKeep = idxKeep[idxKeep > round(0.5*Nrun)]

#idxKeep = 1:500

plotSPTM(ResGibbsM, "hyperparameter", keepRun = idxKeep)
plotSPTM(ResGibbsM, "residuals", keepRun = idxKeep)
plotSPTM(ResGibbsM, "covariates", keepRun = idxKeep)
plotSPTM(ResGibbsM, "cluster", keepRun = idxKeep)



predData = data_plot
predData$yP =  predictSPTM(ResGibbsM, keepRun = idxKeep)$Yp
predData$CI025 = apply(ResGibbsM$GibbsOut$yHat[idxKeep,,idxBatch],2,quantile,0.025)
predData$CI975 = apply(ResGibbsM$GibbsOut$yHat[idxKeep,,idxBatch],2,quantile,0.975)
predData$yPp =   apply(ResGibbsM$GibbsOut$yHat[idxKeep,,idxBatch],2,quantile, 0.5)
predData$yPp =   apply(ResGibbsM$GibbsOut$yHat[idxKeep,,idxBatch],2,mean)

predData = merge(predData, Cldf, by  ="ID", all.x = TRUE)

ggplot(data = predData[predData$Cluster %in% unique(Cldf$Cluster),], aes(x = date, y = obs, col = ID)) + geom_point() + geom_line() + 
  geom_line(aes(x = date, y = yP, group = ID), col = "black") +
  geom_line(aes(x = date, y = yPp, group = ID), col = "red") +
  facet_wrap(Cluster~ID) + 
  geom_ribbon(aes(ymin = CI025, ymax = CI975), fill = "grey70", col = "grey70", alpha = 0.2)

