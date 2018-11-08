library(LearnBayes)
library(mvtnorm)
library(RColorBrewer)
library(bayesm)
library(MASS)
library(compositions)
library(splines)
library(wrapr)
library(slam)
library("foreach")
library("doParallel")
library(Rcpp)
library(tmvtnorm)
library(deldir)
library(tictoc)
library(EDISON)
library(SciencesPo)
load("data/Rdata/dataNitrogenPost2008.Rdata")
source("functions.R")
sourceCpp("src/cpp_functions.cpp")
# Log or not log
#[-c(1,2,8,9,10,11,13)]
subWaterQ.df = droplevels(subset(subWaterQ.df, 
                                 ID %in% levels(subWaterQ.df$ID)[-3]))
subWaterQ.df.TEST = droplevels(subset(subWaterQ.df, date > "2012-01-01"))
subWaterQ.df = droplevels(subset(subWaterQ.df, date <= "2012-01-01"))
subWaterQ.df = na.omit(subWaterQ.df)
subWaterQ.df$obs = log(subWaterQ.df$obs)

subCoords.df = droplevels(subset(subCoords.df, ID %in% levels(subWaterQ.df$ID)))
subLandUse.df$ID = rownames(subLandUse.df)
subLandUse.df = droplevels(subset(subLandUse.df, ID %in% levels(subWaterQ.df$ID), 
                           select = -ID))

ggplot(data = subWaterQ.df, aes(x = date, y = obs, col = ID)) + geom_point() + geom_line() + 
  facet_wrap(~ID, scales = "free") +
  geom_point(data = subWaterQ.df.TEST, aes(x = date, y = log(obs)), col = "black", size=0.5)

# Move on

df = SPTMData(df.obs = subWaterQ.df, tempBasis = "bs", tempPeriod = c("%m", "%Y"), nSplines = c(26,5))
subCoords.df$site = subCoords.df$ID
CLR = clr(subLandUse.df)

SpTcov = cbind(CLR[match(df$obs.data$ID,rownames(subLandUse.df)),])
df.sptmod = SPTModel(df.sptm = df, coordinates = subCoords.df, SpTcov = SpTcov)
set.seed(1)

ResGibbs = estimGibbs(df.sptmod, priors = list(beta = list(m0 = 1, s0 = 0.1, dist = "gauss", constrained = F),
                                          alpha = list(m0 = 0, s0 = 2, dist = "gauss", constrained = F),
                                          gamma = list(m0 = 0, s0 = 0.1, dist = "gauss"),
                                          sigma2 = list(a0=2, b0=1, dist = "igamma"),
                                          tau = list(a0 = 5, b0= 1, dist = "igamma"),
                                          phi = list(inf = 0.1, sup = 50, dist = "unif"),
                                          pi = list(alpha0= matrix(1, nrow=nlevels(df$obs.data$ID), ncol=1), #
                                                    dist = "dirichlet"),
                                          tauT = list(a0 = 2, b0= 1, dist = "gamma"),
                                          phiT = list(inf = 0.25, sup = 0.5, dist = "unif"),
                                          rho=list(inf=2, sup = 3, dist="unif")),
                 N.run = 15000, debug =T,tempRE = "notcorr", print.res = T, 
                 nBatch = 3, parallel = F, nCluster = 1)

saveRDS(object = ResGibbs, file = "data/rds/CLR_results_15000runs_3batches.rds")

df.sptmod = SPTModel(df.sptm = df, coordinates = subCoords.df, SpTcov = NULL)
pi0 = subLandUse.df
pi0 = pi0#[,-c(4:6)]/ rowSums(pi0[,-c(4:6)])
set.seed(1)
ResGibbsM = estimGibbs(df.sptmod, priors = list(beta = list(m0 = 1, s0 = 0.1, dist = "gauss", constrained =F),
                                           alpha = list(m0 = 0, s0 = 2, dist = "gauss", constrained = F),
                                           gamma = list(m0 = 0, s0 = 0.1, dist = "gauss"),
                                           sigma2 = list(a0=2, b0=1, dist = "igamma"),
                                           tau = list(a0 = 5, b0= 1, dist = "igamma"),
                                           phi = list(inf = 0.1, sup = 50, dist = "unif"),
                                           tauT = list(a0 = 2, b0= 1, dist = "gamma"),
                                           phiT = list(inf = 0.1, sup = 0.8, dist = "unif"),
                                           pi = list(alpha0 = pi0, dist = "dirichlet"),
                                           rho=list(inf=2, sup = 3, dist="unif")),
                  N.run = 15000, debug = F, tempRE = "notcorr", model = "simpleMixture",
                  print.res = T, nBatch = 3,
                  parallel = F, nCluster = 7)

saveRDS(object = ResGibbsM, file = "data/rds/Mixture_results_15000runs_3batches.rds")

df.sptmod = SPTModel(df.sptm = df, coordinates = subCoords.df, SpTcov = NULL)
pi0 = subLandUse.df
pi0 = pi0/ rowSums(pi0)
set.seed(1)
ResGibbsMS = estimGibbs(df.sptmod, priors = list(beta = list(m0 = 1, s0 = 0.1, dist = "gauss", constrained =F),
                                                 alpha = list(m0 = 0, s0 = 2, dist = "gauss", constrained = F),# converged with F
                                                 gamma = list(m0 = 0, s0 = 0.1, dist = "gauss"),
                                                 sigma2 = list(a0=2, b0=1, dist = "gamma"),
                                                 tau = list(a0 = 5, b0= 1, dist = "gamma"),
                                                 phi = list(inf = 0.1, sup = 50, dist = "unif"),
                                                 tauT = list(a0 = 2, b0= 1, dist = "gamma"),
                                                 phiT = list(inf = 0.1, sup = 0.8, dist = "unif"),
                                                 pi = list(alpha0 = pi0, dist = "dirichlet"),
                                                rho=list(inf=2, sup = 3, dist="unif")),
                       N.run = 15000, debug = F, model = "spatialMixture",tempRE = "notcorr", 
                       print.res = T, nBatch = 3,
                       parallel = F, nCluster = 7)

saveRDS(object = ResGibbsMS, file = "data/rds/SpatialMixture_results_15000runs_3batches.rds")
