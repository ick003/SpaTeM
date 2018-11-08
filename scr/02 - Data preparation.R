# Preping the data for estimation
load("./data/Rdata/GoyderData.Rdata")
library(reshape2)
# Water measurements

subWaterQ.df = subset(WaterQ.df, param == levels(WaterQ.df$param)[6], select = c(site, date, value))

subWaterQ.df = subset(subWaterQ.df, date < "2010-01-01", select = c(site, date, value))
subWaterQ.df = subset(subWaterQ.df, date > "2008-01-01", select = c(site, date, value))

names(subWaterQ.df) <- c("ID" ,"date", "obs")
subWaterQ.df$ID <- as.factor(as.character(subWaterQ.df$ID))
subWaterQ.df$obs[subWaterQ.df$obs %in% c(0, 0.003, 0.005)] <- NA
#subWaterQ.df$logObs <- log(subWaterQ.df$obs)

ggplot(data = subWaterQ.df, aes(x=date, y = log(obs), col=ID))+geom_point()+geom_line()

subWaterQ.df <- droplevels(subWaterQ.df[-which(subWaterQ.df$ID == "A5051005"),])

subWaterQ.df$ID = as.factor(as.character(subWaterQ.df$ID))

# Land use information

SummaryLU = dcast(LandUse.df[c("DES_SEC_V6","SiteID","Shape_Area")], SiteID ~DES_SEC_V6, sum)

subLandUse.df = subset(SummaryLU, SiteID %in% unique(subWaterQ.df$ID))

#subLandUse.df = rbind(subLandUse.df,1)
#subLandUse.df$SiteID[18] <- as.character(unique(subWaterQ.df$ID)[18])

# One station is missing LU information


siteNames = subLandUse.df$SiteID
rownames(subLandUse.df) <- subLandUse.df$SiteID
subLandUse.df$SiteID <- NULL

subLandUse.df = diag(rowSums(subLandUse.df)^-1) %*% as.matrix(subLandUse.df) 

keep_lu = names(which(apply(subLandUse.df,2,max) > .4))

# Merging GMP and OMU as they seem to behave identically
subLandUse.df = as.data.frame(subLandUse.df)
rownames(subLandUse.df) <- siteNames
subLandUse.df = subset(subLandUse.df, select = keep_lu)
subLandUse.df$`Grazing modified pastures` = subLandUse.df$`Grazing modified pastures` + subLandUse.df$`Other minimal uses`
subLandUse.df$`Services, Transport and Comm` = subLandUse.df$`Services` + subLandUse.df$`Transport and communication`
subLandUse.df$`Other minimal uses` <- NULL
subLandUse.df$`Transport and communication` <- NULL
subLandUse.df$`Services` = NULL

subCoords.df = unique(subset(Coords.df, site %in% rownames(subLandUse.df)))
names(subCoords.df)[1] <- c("ID")

subCoords.df$ID <- as.factor(as.character(subCoords.df$ID))

save(subCoords.df, subLandUse.df, subWaterQ.df, file = "data/Rdata/dataNitrogenPre2010.Rdata")
