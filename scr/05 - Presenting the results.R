#### Tabulated results ####
nRuns = dim(ResGibbs$GibbsOut$theta$alphaH)[1]
nChains = dim(ResGibbs$GibbsOut$theta$alphaH)[4]

PointEst = t(round(rbind(apply(ResGibbs$GibbsOut$theta$gammaH[round(0.5*nRuns):nRuns,,], 2, mean),
                       apply(ResGibbs$GibbsOut$theta$gammaH[round(0.5*nRuns):nRuns,,], 2, quantile, c(0.025, 0.975))),3))
rownames(PointEst) <- c(colnames(CLR))
colnames(PointEst)[1] <- "Post. mean"


keepRun = 3500:5000
basisSplines = do.call(cbind,ResGibbsM$GibbsOut$basis$splines)
basis = ResGibbsM$GibbsOut$basis
nSplines = ncol(basisSplines)
nBasis = 2;nComp = 7

PointEstM = NULL
for(k in 1:nComp){
PointEsttmp = 0
  for(i in 1:nBasis){
    t0 = ResGibbsM$GibbsOut$y$date
    if(basis$tempPeriod[i] == "%m"){
      t0 = as.Date(paste("2014-",format(t0,"%m-%d"), sep = ""))
    }
    xPred = seq(min(t0), max(t0), "week")
    xBasis = predict(basis$splines[[i]],xPred)
    yPred = xBasis %*% t(ResGibbsM$GibbsOut$theta$alphaH[keepRun,basis$idx[[i]],k,1])
    PointEsttmp = PointEsttmp + mean(yPred)
  }
PointEstM = c(PointEstM, PointEsttmp)
}

PointEst = cbind(PointEst, PointEstM)

basisSplines = do.call(cbind,ResGibbsMS$GibbsOut$basis$splines)
basis = ResGibbsMS$GibbsOut$basis
nSplines = ncol(basisSplines)
nBasis = 2;nComp = 7

PointEstMS = NULL
for(k in 1:nComp){
  PointEsttmp = 0
  for(i in 1:nBasis){
    t0 = ResGibbsMS$GibbsOut$y$date
    if(basis$tempPeriod[i] == "%m"){
      t0 = as.Date(paste("2014-",format(t0,"%m-%d"), sep = ""))
    }
    xPred = seq(min(t0), max(t0), "week")
    xBasis = predict(basis$splines[[i]],xPred)
    
    yPred = xBasis %*% t(ResGibbsMS$GibbsOut$theta$alphaH[keepRun,basis$idx[[i]],k,1])
    
    PointEsttmp = PointEsttmp + mean(yPred)
  }
  PointEstMS = c(PointEstMS, PointEsttmp)
}

PointEst = cbind(PointEst, PointEstMS)

round(apply(ResGibbsM$GibbsOut$theta$zH[3500:5000,,,1],2:3, mean),2)
round(apply(ResGibbsMS$GibbsOut$theta$zH[3500:5000,,,1],2:3, mean),2)

#### Plotted results ####

library(RColorBrewer)

plotSPTM(ResGibbs, "tempBasis")
plotSPTM(ResGibbs, "alpha")
plotSPTM(ResGibbs, "beta")
plotSPTM(ResGibbs, "spatialBasis")
plotSPTM(ResGibbs, "hyperparameter", keepRun = 1500:2500)
plotSPTM(ResGibbs, "residuals")
plotSPTM(ResGibbs, "covariates")

predData = df$obs.data
predData$yP =  colMeans(ResGibbs$GibbsOut$yHat[-1,,1])
predData$CI025 = apply(ResGibbs$GibbsOut$yHat[-1,,1],2,quantile,0.025)
predData$CI975 = apply(ResGibbs$GibbsOut$yHat[-1,,1],2,quantile,0.975)

ggplot(data = predData, aes(x = date, y = obs, col = ID)) + geom_point(size=0.5) + #geom_line(size = 0.5) + 
  geom_segment(aes(x = date, xend = date, y = yP, yend = obs))+
  geom_line(aes(x = date, y = yP, group = ID), col = "black") + facet_wrap(~ID, scales = "free_x") + 
  geom_ribbon(aes(ymin = CI025, ymax = CI975), fill = "grey70", col = "grey70", alpha = 0.2)

SSE_S = sqrt(sum((predData$obs - predData$yP)^2))

# Results Mixture

nRun = nrow(ResGibbsM$GibbsOut$theta$rhoH)
Cldf = data.frame(ID = rownames(pi0),Cluster = apply(t(apply(ResGibbsM$GibbsOut$theta$zH[round(3*nRun/5):nRun,,,1],c(2,3),mean)),1,which.max))
idxKeep = which(colSums(apply(ResGibbsM$GibbsOut$theta$zH[,,,1], 1, function(x) apply(x, 2,which.max) - Cldf$Cluster))==0)
idxKeep = idxKeep[idxKeep > 3.*nRun/5]

plotSPTM(ResGibbsM, "tempBasis", keepRun = idxKeep)
plotSPTM(ResGibbsM, "spatialBasis")
plotSPTM(ResGibbsM, "alpha")
plotSPTM(ResGibbsM, "beta")
plotSPTM(ResGibbsM, "hyperparameter", keepRun = idxKeep)
plotSPTM(ResGibbsM, "residuals")
plotSPTM(ResGibbsM, "covariates", keepRun = idxKeep)
plotSPTM(ResGibbsM, "cluster", keepRun =idxKeep)
plotSPTM(ResGibbsM, "map")

# Identify the MCMC samples for which the component allocation corresponds to the most likely one:

predData = df$obs.data
predData$yP =  colMeans(ResGibbsM$GibbsOut$yHat[idxKeep,,1])
predData$CI025 = apply(ResGibbsM$GibbsOut$yHat[idxKeep,,1],2,quantile,0.025)
predData$CI975 = apply(ResGibbsM$GibbsOut$yHat[idxKeep,,1],2,quantile,0.975)
predPost = predictSPTM(ResGibbsM, keepRun = idxKeep, posterior = F)
predData$yPp =  predPost$YpwN
predData$CI025p = predPost$CI025
predData$CI975p = predPost$CI975
predData = merge(predData, Cldf, by  ="ID", all.x = TRUE)

SSE_M = sqrt(sum((predData$obs - predData$yP)^2))

ggplot(data = predData[predData$Cluster %in% 7,], aes(x = date, y = obs, col = ID)) + geom_point(size=0.5) +
  geom_segment(aes(x = date, xend = date, y = yP, yend = obs))+
  geom_line(aes(x = date, y = yP, group = ID), col = "black") + 
  geom_line(aes(x = date, y = yPp, group = ID), col = "blue") + 
  facet_wrap(Cluster~ID, scales = "free") + 
  geom_ribbon(aes(ymin = CI025, ymax = CI975), fill = "grey70", col = "grey70", alpha = 0.2)+
  geom_ribbon(aes(ymin = CI025p, ymax = CI975p), fill = "cyan", col = NA, alpha = 0.2)


nid = "A5040508"
tPred = seq(min(df$obs.data$date), max(df$obs.data$date), by = 'week')
newData = list(ID = rep(as.factor(nid),length(tPred)), 
              time = tPred, 
              loc = subCoords.df[which(subCoords.df$site == nid), 2:3], 
              cov = 1)

rownames(newData$loc) <- nid

ttPred = predictSPTM(ResGibbsM, keepRun = idxKeep, newdata = newData)
# Results Spatial Mixture

testPred = data.frame(yP = ttPred$Yp, 
                      yPp = ttPred$YpwN,
                      date = newData$time, 
                      CI025 = ttPred$CI025, 
                      CI975 = ttPred$CI975)

ggplot(data = predData[predData$Cluster %in% 7,], aes(x = date, y = obs, col = ID)) + 
  geom_segment(aes(x = date, xend = date, y = yP, yend = obs))+
  geom_line(aes(x = date, y = yP, group = ID), col = "black") + 
  geom_line(data = testPred, aes(x = date, y = yPp), col = "blue") +
  facet_wrap(Cluster~ID, scales = "free")+
  geom_ribbon(data = testPred, aes(x = date, y = yP, ymin = CI025, ymax = CI975), fill = "cyan", col = NA,alpha = 0.2) +
  geom_point(size=0.5)
  

nRun = nrow(ResGibbsM$GibbsOut$theta$rhoH)
plotSPTM(ResGibbsMS, "tempBasis")
plotSPTM(ResGibbsMS, "spatialBasis")
plotSPTM(ResGibbsMS, "hyperparameter", keepRun =round(4*nRun/5):nRun)
plotSPTM(ResGibbsMS, "prediction")
plotSPTM(ResGibbsMS, "residuals")
plotSPTM(ResGibbsMS, "covariates")
plotSPTM(ResGibbsMS, "cluster", keepRun =round(4*nRun/5):nRun)
plotSPTM(ResGibbsMS, "map")

predData = df$obs.data
predData$yP =  colMeans(ResGibbsMS$GibbsOut$yHat[-1,,1])
predData$CI025 = apply(ResGibbsMS$GibbsOut$yHat[-1,,1],2,quantile,0.025)
predData$CI975 = apply(ResGibbsMS$GibbsOut$yHat[-1,,1],2,quantile,0.975)
predData$yPp =  colMeans(predictSPTM(ResGibbsM, keepRun = round(4*nRun/5):nRun)$YpwN)
Cldf = data.frame(ID = rownames(pi0),Cluster = apply(t(apply(ResGibbsMS$GibbsOut$theta$zH[round(4*nRun/5):nRun,,,1],c(2,3),mean)),1,which.max))

predData = merge(predData, Cldf, by  ="ID", all.x = TRUE)

SSE_MS = sqrt(sum((predData$obs - predData$yP)^2))

ggplot(data = predData[predData$Cluster %in% c(1:9),], aes(x = date, y = obs, col = ID)) + geom_point(size = 0.5) + 
#  geom_line() + 
  geom_line(aes(x = date, y = yP, group = ID), col = "black") + 
  geom_segment(aes(x = date, xend = date, y = yP, yend = obs))+
  facet_wrap(Cluster~ID, scales = "free") + 
  geom_ribbon(aes(ymin = CI025, ymax = CI975), fill = "grey70", col = "grey70", alpha = 0.2)






