# Using Bayes Factor to select the best model 
# (Following Robert's book [Robert 2012])
library(ggplot2)
library(gridExtra)
library(reshape2)
library(plyr)
library(bayesm)
source("functions.R")

ResGibbs = readRDS("data/rds/CLR_results_1500runs_2batches_constrainedAlpha.rds")
ResGibbsM = readRDS("data/rds/Mixture_results_7500runs_2batches_constrainedAlpha.rds")
ResGibbsMS = readRDS("data/rds/SpatialMixture_results_7500runs_2batches.rds")

nSample = 7500
keepRunS = 1000:1500
keepRunM = 3500:5000


BF_MoS = BF_MSoS =BF_MoMS = BFt = NULL
for(rep in 1:1000){
  print(rep)
mLSimple = margLikelihood(SPTMresobj = ResGibbs, nSample = nSample, keepRun = keepRunS, c = 1)
mLMixture = margLikelihood(SPTMresobj = ResGibbsM, nSample = nSample, keepRun = keepRunM, c = 1)
mLSpatialMixture = margLikelihood(SPTMresobj = ResGibbsMS, nSample = nSample, keepRun = keepRunM, c = 1)

BFt = rbind(BFt, cbind(log(mean(exp(mLSpatialMixture$dPost - max(mLSpatialMixture$dPost))))+ max(mLSpatialMixture$dPost), 
                       log(mean(exp(mLMixture$dPost - max(mLMixture$dPost))))+ max(mLMixture$dPost), 
                       log(mean(exp(mLSimple$dPost - max(mLSimple$dPost))))+ max(mLSimple$dPost)))
BF_MoS = c(BF_MoS, mean(exp(mLMixture$dPost - max(mLMixture$dPost))) / mean(exp(mLSimple$dPost - max(mLSimple$dPost))) * exp(max(mLMixture$dPost)-max(mLSimple$dPost)))
BF_MSoS = c(BF_MSoS, mean(exp(mLSpatialMixture$dPost - max(mLSpatialMixture$dPost))) / mean(exp(mLSimple$dPost - max(mLSimple$dPost))) * exp(max(mLSpatialMixture$dPost)-max(mLSimple$dPost)))
BF_MoMS = c(BF_MoMS, mean(exp(mLMixture$dPost - max(mLMixture$dPost))) / mean(exp(mLSpatialMixture$dPost - max(mLSpatialMixture$dPost))) * exp(max(mLMixture$dPost)-max(mLSpatialMixture$dPost)))
}

saveRDS(BFt, file = "data/rds/BayesFactor_7500_262.rds")

BFr_SMoM = mean(exp(BFt[,1] - max(BFt[,1]))) / mean(exp(BFt[,2] - max(BFt[,2]))) * exp(max(BFt[,1])-max(BFt[,2]))
BFr_SMoS = mean(exp(BFt[,1] - max(BFt[,1]))) / mean(exp(BFt[,3] - max(BFt[,3]))) * exp(max(BFt[,1])-max(BFt[,3]))
BFr_SoM = mean(exp(BFt[,3] - max(BFt[,3]))) / mean(exp(BFt[,2] - max(BFt[,2]))) * exp(max(BFt[,3])-max(BFt[,2]))

postD_SMoM = mean(exp(ResGibbsMS$GibbsOut$logPostDist$Post[keepRunM,,1] - max(ResGibbsMS$GibbsOut$logPostDist$Post[keepRunM,,1], na.rm=T)), na.rm=T) / 
  mean(exp(ResGibbsM$GibbsOut$logPostDist$Post[keepRunM ,,1] - max(ResGibbsM$GibbsOut$logPostDist$Post[keepRunM ,,1], na.rm=T)), na.rm=T) * 
  exp(max(ResGibbsMS$GibbsOut$logPostDist$Post[keepRunM ,,1], na.rm= T )-max(ResGibbsM$GibbsOut$logPostDist$Post[keepRunM ,,1], na.rm=T))

postD_SMoS = mean(exp(ResGibbsMS$GibbsOut$logPostDist$Post[keepRunM,,1] - max(ResGibbsMS$GibbsOut$logPostDist$Post[keepRunM,,1], na.rm=T)), na.rm=T) / 
  mean(exp(ResGibbs$GibbsOut$logPostDist$Post[keepRunS ,,1] - max(ResGibbs$GibbsOut$logPostDist$Post[keepRunS ,,1], na.rm=T)), na.rm=T) * 
  exp(max(ResGibbsMS$GibbsOut$logPostDist$Post[keepRunM ,,1], na.rm= T )-max(ResGibbs$GibbsOut$logPostDist$Post[keepRunS ,,1], na.rm=T))

postD_MoS = mean(exp(ResGibbsM$GibbsOut$logPostDist$Post[keepRunM,,1] - max(ResGibbsM$GibbsOut$logPostDist$Post[keepRunM,,1], na.rm=T)), na.rm=T) / 
  mean(exp(ResGibbs$GibbsOut$logPostDist$Post[keepRunS ,,1] - max(ResGibbs$GibbsOut$logPostDist$Post[keepRunS ,,1], na.rm=T)), na.rm=T) * 
  exp(max(ResGibbsM$GibbsOut$logPostDist$Post[keepRunM ,,1], na.rm= T )-max(ResGibbs$GibbsOut$logPostDist$Post[keepRunS ,,1], na.rm=T))

rangeY = quantile(c(log(c(BF_MoS[is.finite(BF_MoS)]^-1, BF_MSoS[is.finite(BF_MSoS)], BF_MoMS[is.finite(BF_MoMS)]))), c(0.001,0.999))

par(mfrow=c(1,3), mar=c(5,4,4,1))
boxplot(log(BF_MoS), main = "M vs S");abline(h=0);abline(h = -log(BFr_SoM), col = "blue");abline(h = log(postD_MoS), lwd = 2, col = "green")
boxplot(log(BF_MSoS), main = "SM vs S");abline(h=0);abline(h = log(BFr_SMoS), col = "blue");abline(h = log(postD_SMoS), lwd = 2, col = "green")
boxplot(-log(BF_MoMS), main = "SM vs M");abline(h=0);abline(h = log(BFr_SMoM), col = "blue");abline(h = log(postD_SMoM), lwd = 2, col = "green")


