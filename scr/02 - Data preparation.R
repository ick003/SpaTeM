# Preping the data for estimation
load("./data/Rdata/GoyderData.Rdata")
library(tidyverse)
# Water measurements

parameterName = "Nitrate + Nitrite as Nitrogen (mg/L)"

subWaterQ.df.TEST <- WaterQ.df %>% filter(param == parameterName & date > "2011-12-31") %>% 
  filter(site != "A5051005") %>%
  select(ID = site, date, obs = value) %>% mutate(ID = factor(ID)) %>% na.omit()

subWaterQ.df <- WaterQ.df %>% filter(param == parameterName & date > "2000-12-31" & date < "2012-01-01") %>% 
  filter(site != "A5051005") %>%
  select(ID = site, date, obs = value) %>% mutate(ID = factor(ID)) %>% na.omit()

 ggplot(data = subWaterQ.df , aes(x=date, y = log(obs), group = ID, col = ID))+
   geom_point(size = 1.5)+geom_line(size = 0.01) + theme_bw()

# Land use information: select second subdivision, sum the surface, keep only sites with nitrogen measures and remove landuse with no surface

subLandUse.df <- LandUse.df %>% select(SiteID, DES_SEC_V6, Shape_Area) %>% group_by(SiteID, DES_SEC_V6) %>% 
  summarise(Shape_Area = sum(Shape_Area)) %>% mutate(Shape_Area = round(Shape_Area / sum(Shape_Area), 3)) %>%  
  filter(Shape_Area > 0.05) %>%
  spread(DES_SEC_V6, Shape_Area) %>% filter(SiteID %in% subWaterQ.df$ID) %>% 
  select_if(function(x){!all(is.na(x))}) %>% 
  replace(is.na(.), 0) %>% 
  mutate(`Grazing modified pastures` = `Grazing modified pastures` + `Other minimal uses`) %>% 
  mutate(`Services, Transport and Comm` = `Services` + `Transport and communication`) %>%
  select(-c(`Services`, `Other minimal uses`, `Transport and communication`)) %>%
  column_to_rownames('SiteID') 
  

names(Coords.df)[4] <- "Else"
subCoords.df <- Coords.df %>% filter(site %in% rownames(subLandUse.df)) %>% select(site, longitude, latitude) %>% 
  distinct() %>% droplevels() 

save(subCoords.df, subLandUse.df, subWaterQ.df, subWaterQ.df.TEST, file = "data/Rdata/dataNitrogenPost2008.Rdata")

save(subCoords.df, subLandUse.df, subWaterQ.df, subWaterQ.df.TEST, file = "data/Rdata/dataNitrogenPost2000.Rdata")
