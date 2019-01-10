# Data section: Map and overview

# Figure 2 - Mount Lofty range and Australia

library(rgdal)
library(tidyverse)
library(sf)
library(ggpubr)
library(gridExtra)
library(rgeos)

theme_set(theme_bw())
theme_opts <- theme(
  text = element_text(color = "#22211d"), 
  plot.title = element_text(size= 22, hjust=0.5, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm")),
  plot.subtitle = element_text(size= 13, hjust=0.5, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm")),
  legend.position = c(0.2, 0.26),
  panel.border = element_blank()
)

#### Map of measurement locations ####

world_GAUL1 = readOGR("./data/shapefiles/ne_50m_admin_1_states_provinces/",    layer="ne_50m_admin_1_states_provinces") 

australia <- st_as_sf(world_GAUL1) %>% filter(admin == "Australia")  %>% mutate_if(is.factor, as.character) 
australia <- cbind(australia, st_coordinates(st_centroid(australia)))

# Correct ACT y cordinates for map readability
australia$Y[australia$name == "Australian Capital Territory"] = -36.0

load("./data/Rdata/GoyderData.Rdata")
proj4string(CatchmentMap) <- proj4string(world_GAUL1)

catchmentBox <- as(raster::extent(CatchmentMap@bbox), "SpatialPolygons")
proj4string(catchmentBox)<- proj4string(world_GAUL1)

# Streams

stream = readOGR("./data/shapefiles/stream_australia/",    layer="au_riv_30s")

stream <- spTransform(stream, CRS(proj4string(CatchmentMap)))
stream_sa <- stream[CatchmentMap,]
river_mountlofty <- st_as_sf(stream_sa)

pAustralia <- ggplot(data = australia) + 
 geom_sf(size = 0.15) + 
  geom_text(aes(x=X, y=Y, label = name), size=2.5)+
  annotate("rect", xmin = catchmentBox@bbox[1,1], xmax = catchmentBox@bbox[1,2], 
           ymin = catchmentBox@bbox[2,1], ymax = catchmentBox@bbox[2,2],alpha = .185, fill = "red", color="red") +
# coord_sf(ylim = c(-44,-9), expand = FALSE) +
  coord_sf(datum= NA)+
  theme_opts + xlab("") + ylab("")

MountLofty <- st_as_sf(CatchmentMap)
sastate_poly <- readOGR("./data/shapefiles/sastatepolygon/", layer="SA_STATE_POLYGON_shp")
sastate <- st_as_sf(sastate_poly)
Stations <- WaterQ.df %>% group_by(site, Nitrogen = param == "Nitrate + Nitrite as Nitrogen (mg/L)") %>% 
  summarise() %>% group_by(site, Nitrogen = sum(Nitrogen)) %>% summarise() %>% left_join(Coords.df[,-4]) %>% na.omit()
Stations$Nitrogen <- as.factor(Stations$Nitrogen)
StationsName = Stations %>% filter(Nitrogen == 1) %>% distinct()

levels(Stations$Nitrogen) <- list("Station w Nitrogen"="1", "Station w/o Nitrogen"="0")

StationsName$latitude[StationsName$site == "A5040508"] = -34.8489
StationsName$latitude[StationsName$site == "A5030509"] = -35.04613
StationsName$latitude[StationsName$site == "A5031006"] = StationsName$latitude[StationsName$site == "A5031006"] - 0.015
StationsName$longitude[StationsName$site == "A5031006"] = 138.8011

StationsName$latitude[StationsName$site == "A5031007"] = StationsName$latitude[StationsName$site == "A5031007"] - 0.015
StationsName$longitude[StationsName$site == "A5031007"] = 138.6611

StationsName$latitude[StationsName$site == "A5031008"] = StationsName$latitude[StationsName$site == "A5031008"] - 0.01
StationsName$longitude[StationsName$site == "A5031008"] = 138.8011

StationsName$latitude[StationsName$site == "A5030526"] = StationsName$latitude[StationsName$site == "A5030526"]
StationsName$longitude[StationsName$site == "A5030526"] = 138.6611

pMountLofty <- 
  ggplot(data = MountLofty) + 
  geom_sf(fill = "grey",aes(colour="Subcatchment"),size = 0.2) + 
  geom_sf(data = sastate,fill = "white",alpha = 0.5) + 
  geom_sf(data = river_mountlofty, aes(colour="Streams"), size = 0.2, show.legend = "line") + 
  geom_point(data = Stations, aes(x = longitude, y = latitude,color = Nitrogen), shape = 5)+
  geom_text(data = StationsName, aes(x = longitude, y = latitude+0.015, label = site), size = 2.5)+
  coord_sf(ylim = catchmentBox@bbox[2,], xlim = catchmentBox@bbox[1,]) + 
  scale_colour_manual(name = NULL, values = c("Station w/o Nitrogen" = "darkgreen","Station w Nitrogen"="red", "Subcatchment" = "white", "Streams" = "blue"),
                      guide = guide_legend(override.aes = list(linetype = c("blank","blank", "solid", "solid"), 
                                                               shape = c(5,5, NA, NA)))) + 
#    scale_fill_manual(name = NULL, values = c("Subcatchment" = "grey")) +
  ggsn::scalebar(MountLofty, dist = 10, st.size=3,
                 st.dist = 0.0006, height=0.0004, dd2km = TRUE, model = 'WGS84')+ xlab("") + ylab("")

legend <- get_legend(pMountLofty)

png(file = "./figure/Map.png", res = 300, units = 'in',
    width = 9, height = 12)
grid.arrange(
  grobs = list(pMountLofty + theme(legend.position="none") , pAustralia, legend),
  widths = c(1, 1),
  layout_matrix = rbind(c(2, 3),
                        c(1, 1),
                        c(1,1))
)
dev.off()

#### Table of number of measures per site ####
library(knitr)
library(kableExtra)

NitrateDF <- WaterQ.df %>% group_by(site, date, Nitrogen = param == "Nitrate + Nitrite as Nitrogen (mg/L)") %>% 
  select(site, date , Nitrogen) %>% filter(Nitrogen == TRUE, date >= "2008-01-01")

kable(table(NitrateDF$site, format(NitrateDF$date, "%Y")), format = "latex", booktabs = T,
      caption = "Number of Nitrate observations each year (post $2008$) at each sampling site.")

LandUseDF <- LandUse.df %>% select(SiteID, DES_SEC_V6, Shape_Area) %>% group_by(SiteID, DES_SEC_V6) %>% 
  summarise(Shape_Area = sum(Shape_Area)) %>% mutate(Shape_Area = round(Shape_Area / sum(Shape_Area), 3)) %>%  
  filter(Shape_Area > 0.05) %>%
  spread(DES_SEC_V6, Shape_Area) %>% filter(SiteID %in% NitrateDF$site) %>% 
  select_if(function(x){!all(is.na(x))})

kable(LandUseDF,format = "latex", booktabs = T, caption = "Compositional landuse information at each sampling site.") %>%
  row_spec(0, angle = 45)
