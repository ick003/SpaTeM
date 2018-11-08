# Using posterior predictive distribution on out-of-sample data to test the models
# (Following Gelman's book  [Gelman et al 2013]) 
# ANZECC guidelines for freshwater aquatic ecosystems h =  0.1 mg / L for public water system
h = 0.1

ResGibbs = readRDS("data/rds/CLR_results_1500runs_2batches_constrainedAlpha.rds")
ResGibbsM = readRDS("data/rds/Mixture_results_7500runs_2batches_constrainedAlpha.rds")
ResGibbsMS = readRDS("data/rds/SpatialMixture_results_7500runs_2batches.rds")

keepRunS = 1000:1500
keepRunM = 3500:5000

keepRunMt = sort(ResGibbsM$GibbsOut$logPostDist$Post[keepRunM,1,1], index.return=T, decreasing = T)$ix[1:500]
keepRunMSt = sort(ResGibbsMS$GibbsOut$logPostDist$Post[keepRunM,1,1], index.return=T, decreasing = T)$ix[1:500]

plot(keepRunM,ResGibbsM$GibbsOut$logPostDist$Post[keepRunM,1,1])
points(keepRunM[keepRunMt], ResGibbsM$GibbsOut$logPostDist$Post[keepRunM,1,1][keepRunMt], col="red")


library(reshape2)
library(ggplot2)
library(plyr)

#### CLR model checking ####

yPred = predictSPTM(ResGibbs, posterior = T, keepRun = keepRunS)
yRep  = as.data.frame(t(yPred$YpwN))

# Test Statistic 1 : #times above threshold h, per site 

yRep$site = ResGibbs$GibbsOut$y$ID
yR = melt(yRep)

yRT = ddply(yR, .(site, variable), summarise, prop = sum(value > log(h)))
y0s <- ddply(ResGibbs$GibbsOut$y, .(ID), summarise, prop = sum(obs > log(h)))

# Test Statistic 2 : max difference to threshold per site 

yRT = merge(yRT, ddply(yR, .(site, variable), summarise, amp = max(value) - log(h)))
y0s <- merge(y0s,ddply(ResGibbs$GibbsOut$y, .(ID), summarise, amp = max(obs) - log(h)))

yRT$model = 'CLR'

#### Mixture model checking ####

yPredM = predictSPTM(ResGibbsM, posterior = T, keepRun = keepRunM[keepRunMt])
yRepM  = as.data.frame(t(yPredM$YpwN))

# Test Statistic 1 : #times above threshold h, per site 

yRepM$site =ResGibbsM$GibbsOut$y$ID
yRM = melt(yRepM)

yRTM = ddply(yRM, .(site, variable), summarise, prop = sum(value > log(h)))
y0sM <- ddply(df$obs.data, .(ID), summarise, prop = sum(obs > log(h)))

# Test Statistic 2 : max difference to threshold per site 

yRTM = merge(yRTM, ddply(yRM, .(site, variable), summarise, amp = max(value) - log(h)))
y0sM <- merge(y0sM,ddply(df$obs.data, .(ID), summarise, amp = max(obs) - log(h)))

yRTM$model = 'mixture'

#### SPatial Mixture model checking ####

yPredMS = predictSPTM(ResGibbsMS, posterior = T, keepRun = keepRunM[keepRunMSt])
yRepMS  = as.data.frame(t(yPredMS$YpwN))

# Test Statistic 1 : #times above threshold h, per site 

yRepMS$site = ResGibbsMS$GibbsOut$y$ID
yRMS = melt(yRepMS)

yRTMS = ddply(yRMS, .(site, variable), summarise, prop = sum(value > log(h)))
y0sMS <- ddply(df$obs.data, .(ID), summarise, prop = sum(obs > log(h)))

# Test Statistic 2 : max difference to threshold per site 

yRTMS = merge(yRTMS, ddply(yRMS, .(site, variable), summarise, amp = max(value) - log(h)))
y0sMS <- merge(y0sMS,ddply(df$obs.data, .(ID), summarise, amp = max(obs) - log(h)))

yRTMS$model = 'spatial mixture'

#### PLotting ####


yRFinal = rbind(yRT, yRTM, yRTMS)
pdf(file = "../PredLike.pdf", width = 10, height = 6)
ggplot(data = yRFinal, aes(x = site, y = prop, colour = model)) + geom_boxplot()+
               stat_summary(data = y0s, aes(x = ID, y = prop), fun.y = "mean",  colour = "red", size = 0.5, geom = "point")+
               geom_boxplot(data = y0s, aes(x = ID, y = prop), col="red", size = 0.5) +
  #scale_colour_grey(start = 0.2, end = 0.8) +
  #scale_colour_brewer(type = "qual", palette = 1, direction = 1, aesthetics = "colour")+ 
  scale_colour_brewer(palette = "Dark2")+
  theme_bw() + labs(x = "Sites", y = "Number of days above threshold") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()
