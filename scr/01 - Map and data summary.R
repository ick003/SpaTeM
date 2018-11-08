# Data section: Map and overview

library(maptools)
library(RColorBrewer)
library(classInt)
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(rasterVis)
library(PBSmapping)
library(GEOmap)
library(mapmisc)

shape_path <- "./data/shapefiles/"
coast_shapefile <- paste(shape_path, "ne_50m_coastline/ne_50m_coastline.shp", sep="")
ocean_shapefile <- paste(shape_path, "ne_50m_ocean/ne_50m_ocean.shp", sep="")
admin0_shapefile <- paste(shape_path, "ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp", sep="")
admin1_shapefile <- paste(shape_path, "ne_50m_admin_1_states_provinces/ne_50m_admin_1_states_provinces.shp", sep="")
lakes_shapefile <- paste(shape_path, "ne_50m_lakes/ne_50m_lakes.shp", sep="")
bb_shapefile <- paste(shape_path, "ne_50m_graticules_all/ne_50m_wgs84_bounding_box.shp", sep="")
grat30_shapefile <- paste(shape_path, "ne_50m_graticules_all/ne_50m_graticules_30.shp", sep="")
cities_shapefile <- paste(shape_path, "ne_10m_populated_places/ne_10m_populated_places.shp", sep="")
sastate_shapefile <- paste(shape_path, "sastatepolygon/SA_STATE_POLYGON_shp.shp", sep="")
stream_shapefile <- paste(shape_path, "stream_australia/au_riv_30s.shp", sep="")
  
layer <- ogrListLayers(coast_shapefile)
coast_lines <- readOGR(coast_shapefile, layer=layer)
layer <- ogrListLayers(ocean_shapefile)
ocean_poly <- readOGR(ocean_shapefile, layer=layer)
layer <- ogrListLayers(admin0_shapefile)
admin0_poly <- readOGR(admin0_shapefile, layer=layer)
layer <- ogrListLayers(admin1_shapefile)
admin1_poly <- readOGR(admin1_shapefile, layer=layer)
layer <- ogrListLayers(lakes_shapefile)
lakes_poly <- readOGR(lakes_shapefile, layer=layer)
lrglakes_poly <- lakes_poly[as.numeric(as.factor(lakes_poly$scalerank)) <= 2 ,]
layer <- ogrListLayers(grat30_shapefile)
grat30_lines <- readOGR(grat30_shapefile, layer=layer)
layer <- ogrListLayers(cities_shapefile)
cities_points <- readOGR(cities_shapefile, layer=layer)
layer <- ogrListLayers(sastate_shapefile)
sastate_poly <- readOGR(sastate_shapefile, layer=layer)
layer <- ogrListLayers(stream_shapefile)
stream_lines <- readOGR(stream_shapefile, layer=layer)

layer <- ogrListLayers(bb_shapefile)
bb_poly <- readOGR(bb_shapefile, layer=layer)
bb_lines <- as(bb_poly, "SpatialLines")

robin_crs <- CRS("+proj=robin +lon_0=0w")

bb_poly_proj <- spTransform(bb_poly, robin_crs)
coast_lines_proj <- spTransform(coast_lines, robin_crs)
ocean_poly_proj <- spTransform(ocean_poly, robin_crs)
admin0_poly_proj <- spTransform(admin0_poly, robin_crs)
admin1_poly_proj <- spTransform(admin1_poly, robin_crs)
lakes_poly_proj <- spTransform(lakes_poly, robin_crs)
cities_point_proj <- spTransform(cities_points, robin_crs)
grat30_lines_proj <- spTransform(grat30_lines, robin_crs)
lrglakes_poly_proj <- spTransform(lrglakes_poly, robin_crs)

# Loading Goyder data

load("./data/Rdata/GoyderData.Rdata")
proj4string(CatchmentMap) <- proj4string(admin1_poly)
CatchmentMap_proj <- spTransform(CatchmentMap, robin_crs)

catchmentBox <- as(raster::extent(CatchmentMap@bbox), "SpatialPolygons")
proj4string(catchmentBox)<- proj4string(admin1_poly)

sastateBox <- as(raster::extent(sastate_poly@bbox), "SpatialPolygons")
proj4string(sastateBox)<- proj4string(admin1_poly)

catchmentBox_proj <- spTransform(catchmentBox, robin_crs)

samplingSites_points <- SpatialPointsDataFrame(coords = Coords.df[,c("longitude","latitude")], data = Coords.df,
                                                      proj4string =CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

australia = admin1_poly[admin1_poly@data$admin == "Australia",]
australia_proj = admin1_poly_proj[admin1_poly_proj@data$admin == "Australia",]

States = intersect(australia,admin0_poly)
States_proj = intersect(australia_proj,admin0_poly_proj)
Cities = intersect(cities_points, sastateBox)

Streams = intersect(stream_lines,sastateBox)

StreamsMLR = intersect(stream_lines,catchmentBox)
#save(StreamsMLR, file = "data/Rdata/StreamPolyMLR.RData")

pngfile <- "./figures/AustraliaMap.png"
png(file=pngfile,width=400,height=450)
par(mar=c(1,1,1,1))
plot(australia, col="gray95", bor = "grey")
plot(ocean_poly, bor=NA, col=adjustcolor("darkblue",alpha.f = 0.1), add=T)
plot(admin0_poly, col=NA, bor="grey", add=TRUE)
plot(catchmentBox, col=NA, bor = "red", add=T)
plot(lakes_poly, bor="lightblue", add=TRUE)
plot(lrglakes_poly, bor="lightblue", add=TRUE)
plot(grat30_lines, col="black",lwd=0.2, add=TRUE)
set.seed(1)
text(coordinates(States)+runif(2*nrow(States@data),-0.2,0.2),labels=toupper(States$name.1), cex=0.5, col="darkgreen")
box(lty = 1, lwd=2)
dev.off()



#### Example of times series for Notrigen data



subWaterQ.df = subset(WaterQ.df, param == levels(WaterQ.df$param)[6], select = c(site, date, value))
subWaterQ.df = subset(subWaterQ.df, date > "2008-01-01", select = c(site, date, value))
library(xtable)
xtable(table(subWaterQ.df$site,format(subWaterQ.df$date, "%Y")))

subWaterQ.df = subset(subWaterQ.df, date < "2013-01-01", select = c(site, date, value))


#ID_keep = names(which(table(subWaterQ.df$site)>30))
#subWaterQ.df = subset(subWaterQ.df, site %in% ID_keep )
names(subWaterQ.df) <- c("ID" ,"date", "obs")
subWaterQ.df$ID <- as.factor(as.character(subWaterQ.df$ID))
subWaterQ.df$obs[subWaterQ.df$obs==0] <- NA
subWaterQ.df$obs <- log(subWaterQ.df$obs)

ID_plot = c(names(which.max(table(subWaterQ.df$ID))), names(which.min(table(subWaterQ.df$ID))))

### Land use data summary
library(reshape2)
SummaryLU = dcast(LandUse.df[c("DES_SEC_V6","SiteID","Shape_Area")], SiteID ~DES_SEC_V6, sum)

subLandUse.df = subset(SummaryLU, SiteID %in% unique(subWaterQ.df$ID))
siteNames = subLandUse.df$SiteID


### PLot Nitrogen stations #
colStation = rep("black", nrow(samplingSites_points@data))
colStation[which(samplingSites_points@data$site %in% siteNames)] = "red"

pngfile <- "./figures/CatchmentMap.png"
png(file=pngfile,width=400,height=450)
par(mar=c(1,1,1,1))
plot(catchmentBox, bg = adjustcolor("darkblue",alpha.f = 0.1), bor = NA)
plot(sastate_poly, add=T, bor="grey", col="white")
plot(CatchmentMap, col="gray95", bor="pink", add=T)
plot(Streams, col="lightblue", add=T)
plot(Cities, col="darkgreen", add=T, cex=0.5, pch=5)
text(coordinates(Cities)+rep(c(0,0.02),each=nrow(Cities@data)),labels=toupper(Cities$NAME), cex=0.6, col="darkgreen")
plot(samplingSites_points, add=T, cex=0.5, col = colStation)
box(lty = 1, lwd=2)
scalebar(20, below="kilometers",type='bar', lonlat=T, cex=0.85, label = c(0,10,20), 
         xy = samplingSites_points@bbox[cbind(c(1,2),c(1,2))])
axis(1, at = seq(137,140,0.3), parse(text=degreeLabelsEW(seq(137,140,0.3))), padj = -5, tck = 0.01, cex.axis=0.65)
axis(2, at = seq(-36,-34,0.2), parse(text=degreeLabelsNS(seq(-36,-34,.2))), padj = 5, tck = 0.01, cex.axis=0.6)
dev.off()

pngfile <- "./figures/LegendMap.png"
png(file=pngfile,width=800,height=100)
par(mar=c(1,1,1,1))
plot(0,0,bty='n', axes=F, xlim = c(1,1.5))
legend("center", legend = c("Stream","Subcatchment boundary","City","Stations w/o Nitrogen", "Stations w Nitrogen"), 
       bty='n',lty = c(1,1,NA,NA,NA),pch = c(NA,NA,5,3,3),cex=1.2,
       col = c("lightblue","pink","darkgreen","black","red"),lwd=c(2,2,NA,NA,NA),
       ncol=3)
box(lty=1, lwd=2)
dev.off()


CatchmentMapD = CatchmentMap[which(!duplicated(CatchmentMap$NewCatch)),]

#### Two-in-one map ####

colStation = rep("darkgreen", nrow(samplingSites_points@data))
colStation[which(samplingSites_points@data$site %in% siteNames)] = "red"

pngfile <- "./figures/CatchmentMap2.png"
png(file=pngfile,width=400,height=450)
par(mar=c(1,1,1,1))
plot(catchmentBox, bg = adjustcolor("white",alpha.f = 0.1), bor = NA)
plot(sastate_poly, add=T, bor="grey", col=adjustcolor("darkgrey", alpha.f = 0.3))
plot(CatchmentMapD, col=adjustcolor("orange", alpha.f = 0.5), bor="grey", add=T)
plot(StreamsMLR, col="lightblue", add=T)
plot(Cities, col="black", add=T, cex=0.5, pch=5)
text(coordinates(Cities)+rep(c(0,0.02),each=nrow(Cities@data)),labels=tolower(Cities$NAME), cex=0.75, col="black")
plot(samplingSites_points, add=T, cex=0.5, col = colStation)
box(lty = 1, lwd=2)
scalebar(20, below="kilometers",type='line', lonlat=T, cex=0.85, label = c("20 km"), 
         xy = samplingSites_points@bbox[cbind(c(1,2),c(1,2))])
axis(1, at = seq(137,140,0.3), parse(text=degreeLabelsEW(seq(137,140,0.3))), padj = -5, tck = 0.01, cex.axis=0.65)
axis(2, at = seq(-36,-34,0.2), parse(text=degreeLabelsNS(seq(-36,-34,.2))), padj = 5, tck = 0.01, cex.axis=0.6)
dev.off()


pngfile <- "./figures/AustraliaMap2.png"
png(file=pngfile,width=400,height=450)
par(mar=c(1,1,1,1))
plot(australia, col="gray95", bor = "grey")
plot(sastate_poly, add=T, bor="grey", col=adjustcolor("darkgrey", alpha.f = 0.3))
plot(ocean_poly, bor=NA, col=adjustcolor("white",alpha.f = 0.1), add=T)
plot(admin0_poly, col=NA, bor="grey", add=TRUE)
plot(catchmentBox, col=NA, bor = "red", add=T, lwd=3)
#plot(lakes_poly, bor="lightblue", add=TRUE)
#plot(lrglakes_poly, bor="lightblue", add=TRUE)
#plot(grat30_lines, col="black",lwd=0.2, add=TRUE)
set.seed(1)
text(coordinates(States)+runif(2*nrow(States@data),-0.2,0.2),labels=toupper(States$name.1), cex=0.5, col="black")
box(lty = 1, lwd=2)
dev.off()


pngfile <- "./figures/LegendMap2.png"
png(file=pngfile,width=200,height=200)
par(mar=c(1,0.2,1,0.2))
plot(0,0,bty='n', axes=F, xlim = c(1,1.5))
legend("center", legend = c("Stream","Subcatchment boundary","City","Stations w/o Nitrogen", "Stations w Nitrogen"), 
       bty='n',lty = c(1,1,NA,NA,NA),pch = c(NA,NA,5,3,3),cex=1,
       col = c("lightblue","grey","black","darkgreen","red"),lwd=c(2,2,NA,NA,NA),
       ncol=1)
#box(lty=1, lwd=2)
dev.off()

# Figure 2 - Mount Lofty range and Australia

library(data.table)
library(rgdal)
library(ggplot2)
library(rgeos)
library(tidyverse)
library(broom)
library(viridis)
library(maptools)   
library(sf)
library(maps)
library(egg)
library(ggsn)

theme_set(theme_bw())
theme_opts <- theme(
  text = element_text(color = "#22211d"), 
  plot.title = element_text(size= 22, hjust=0.5, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm")),
  plot.subtitle = element_text(size= 13, hjust=0.5, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm")),
  legend.position = c(0.2, 0.26),
  panel.border = element_blank()
)

world_GAUL1 = readOGR("../JABES_paper/data/shapefiles/ne_50m_admin_1_states_provinces/",    layer="ne_50m_admin_1_states_provinces") 

australia <- st_as_sf(world_GAUL1) %>% filter(admin == "Australia")  %>% mutate_if(is.factor, as.character) 
australia <- cbind(australia, st_coordinates(st_centroid(australia)))

load("./data/Rdata/GoyderData.Rdata")
proj4string(CatchmentMap) <- proj4string(world_GAUL1)

catchmentBox <- as(raster::extent(CatchmentMap@bbox), "SpatialPolygons")
proj4string(catchmentBox)<- proj4string(world_GAUL1)

pAustralia <- ggplot(data = australia) + 
 geom_sf() + 
  geom_text(aes(x=X, y=Y, label = name), size=2.5)+
  annotate("rect", xmin = catchmentBox@bbox[1,1], xmax = catchmentBox@bbox[1,2], 
           ymin = catchmentBox@bbox[2,1], ymax = catchmentBox@bbox[2,2],alpha = .185, fill = "red", color="red") +
# coord_sf(ylim = c(-44,-9), expand = FALSE) +
  coord_sf(datum= NA)+
  theme_opts + xlab("") + ylab("")

MountLofty <- st_as_sf(CatchmentMap)
sastate_poly <- readOGR("../JABES_paper/data/shapefiles/sastatepolygon/", layer="SA_STATE_POLYGON_shp")
sastate <- st_as_sf(sastate_poly)
Stations <- WaterQ.df %>% group_by(site, Nitrogen = param == "Nitrate + Nitrite as Nitrogen (mg/L)") %>% 
  summarise() %>% group_by(site, Nitrogen = sum(Nitrogen)) %>% summarise() %>% left_join(Coords.df[,-4]) %>% na.omit()
Stations$Nitrogen <- as.logical(Stations$Nitrogen)
StationsName = Stations %>% filter(Nitrogen == TRUE) %>% distinct()

#pMountLofty <- 
  ggplot(data = MountLofty) + 
  geom_sf(fill = "orange", size = 0.2) + 
  geom_sf(data = sastate, alpha = 0.5) + 
  geom_point(data = Stations, aes(x = longitude, y = latitude, color = Nitrogen), shape = 5)+
  geom_text(data = StationsName, aes(x = longitude, y = latitude+0.025, label = site), size = 2.5)+
  coord_sf(ylim = catchmentBox@bbox[2,], xlim = catchmentBox@bbox[1,]) + 
  scale_colour_manual(values = c("darkgreen", "red")) + theme(legend.position="none") + 
  ggsn::scalebar(MountLofty, dist = 10, st.size=2.5,
                 st.dist = 0.0004, height=0.0003, dd2km = TRUE, model = 'WGS84')

 
grid.arrange(
  grobs = list(pMountLofty, pAustralia),
  widths = c(1, 1, 1.2),
  layout_matrix = rbind(c(1, 1, 2),
                        c(1, 1,NA))
)
