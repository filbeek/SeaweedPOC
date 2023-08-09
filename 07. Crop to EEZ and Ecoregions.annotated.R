######################################################################################
# this script combines the distribution of marine forests but uses the NPP and export rasters.
# it outputs files that measure the area of seaweedforest and the average and sd of export for every EEZ and ocean basin
# the output file is called 'resultsExport' 
# written by K Filbee-Dexter and A. Pessarrodona  for the EUROMARINE project 21.3.22

######################################################################################

# Clear workspace for a faster run -----------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
options("install.lock"=FALSE)

# Libraries -----------------------------------
library(raster)
library(sf)
#devtools::install_github("mdsumner/sfraster")
library(sfraster)
library(sdmpredictors)
library(readxl)
library(tidyverse)
library(ggrepel)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(leaflet)
library(rgdal)
library(maptools)
library(rgeos)
library(meowR)

setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")
rasterOptions(tmpdir="tempDir/")

# Load depth data -----------------------------------
#get depth data, and crop to 0 to 50
#no need to do this!

bathymetryLayer=raster('Datasets/bathy.raster.tif')
#plot(bathy.raster)

###### ...Subset by depth ----
## crop by 50 m value.

bathymetryLayer[bathymetryLayer<=-100]<-NA ## all the other values give them NA


# load export layer and crop it to bathymetry--------

#Think we can remove the next lines
#load export layer and crop it to bathymetry. 
setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")
rasterLayer <- raster("Outputs/raster.percent.export.dist.sinking.2023.11.tif")
rasterLayerF <- raster("Outputs/raster.percent.export.dist.floating.2023.11.tif")


#check layers
#rasterLayer[rasterLayer == 0] <- NA
#rasterLayer[!is.na(rasterLayer)] <- 1

#sinking
rasterLayerHR <- resample(rasterLayer,bathymetryLayer, method="ngb")
rasterLayerHR <- raster::mask(rasterLayerHR,bathymetryLayer)

#floating
rasterLayerHRF <- resample(rasterLayerF,bathymetryLayer, method="ngb")
rasterLayerHRF <- raster::mask(rasterLayerHRF,bathymetryLayer)
#writeRaster(rasterLayerHR,filename="Data/export.tif",format="GTiff",overwrite=T)


#load gC exported layer for floating and sinking species and crop it to bathymetry. 
rasterLayergc.F <- raster("Outputs/raster.gCexported.2023.10.floating.tif")
rasterLayergc.S <- raster("Outputs/raster.gCexported.2023.10.sinking.tif")

rasterLayerGC.S <- resample(rasterLayergc.S,bathymetryLayer, method="ngb")
rasterLayerGC.S <- raster::mask(rasterLayerGC.S,bathymetryLayer)

rasterLayerGC.F <- resample(rasterLayergc.F,bathymetryLayer, method="ngb")
rasterLayerGC.F <- raster::mask(rasterLayerGC.F,bathymetryLayer)

mean(values(rasterLayerGC.F), na.rm = T)
mean(values(rasterLayerGC.S), na.rm = T)

# load shapefiles to get boundaries for ocean basins and the worlds EEZs--------
# no need to run every time
# this can take a while to run
oceanBasins <- shapefile("Global Oceans/goas_v01.shp")
eez <- shapefile("EEZ/EEZGlobal.shp")

sf::sf_use_s2(FALSE)

oceanBasins <- st_as_sf(oceanBasins)
oceanBasins <- st_simplify(oceanBasins, preserveTopology=TRUE, dTolerance = 0.001)
oceanBasins <- st_make_valid(oceanBasins)

eez <- st_as_sf(eez)
eez <- st_simplify(eez, preserveTopology=TRUE, dTolerance = 0.001)
eez <- st_make_valid(eez)

#save these simplified files so you do not have to load for later runs
#st_write(eez, "Outputs/EEZ.csv", layer_options = "GEOMETRY=AS_XY")
#st_write(eez, "Outputs/EEZ.asWKT.csv", layer_options = "GEOMETRY=AS_WKT")
#st_write(oceanBasins, "Outputs/oceanBasins.csv", layer_options = "GEOMETRY=AS_XY")
#st_write(oceanBasins, "Outputs/oceanBasinsas.WKT.csv", layer_options = "GEOMETRY=AS_WKT")

#eez=st_read("Outputs/EEZ.csv", layer_options = "GEOMETRY=AS_XY")
#oceanBasins=st_read("Outputs/oceanBasins.csv")


# -------------------------------------------------#
#### Calculate export for each EEZ area  ####               
#this calculates the average percent export and gC exported for each EEZ WITHOUT dividing by ocean basin 
#it only runs EEZs with seaweed areas. 1:23,25:174,178:186, 191:278
#these calculations are done separately for both floating and sinking species
# -------------------------------------------------#

#floating species 
Subtidal.brown.NPP <- raster("Datasets/Subtidal BrownsPresent.tif") 
Subtidal.brown.NPP.resampled= resample(Subtidal.brown.NPP, # the raster you want to fix the resolution of (and keep the values)
                                       bathymetryLayer) #the other raster
Subtidal.brown.NPP.resampled <- clamp(Subtidal.brown.NPP.resampled, -Inf, 1.684)
max(values(Subtidal.brown.NPP.resampled), na.rm=T)

resultsExportEEZOnly= data.frame()

for( i in c(1:23,25:174,178:186, 191:278) ) {
  
  cat("\014")
  cat("\n")
  cat("# ---------------------------------\n")
  cat("i: ",i,"out of",nrow(eez),"\n")
  cat("# ---------------------------------\n")
  cat("\n")
  
  #i=32
  mask <- eez[i,]
  #floating export
  rasterLayerHRF.i <- crop(rasterLayerHRF,extent(mask))
  cells <- sfraster::cellFromPolygon(rasterLayerHRF.i, mask[,1])
  exportEstimate.F <- mean(rasterLayerHRF.i[unlist(cells)], na.rm=T)
  exportSD.F <- sd(rasterLayerHRF.i[unlist(cells)], na.rm=T)
  #sinking export
  rasterLayerHRS.i <- crop(rasterLayerHR,extent(mask))
  cells <- sfraster::cellFromPolygon(rasterLayerHRS.i, mask[,1])
  exportEstimate.S <- mean(rasterLayerHRS.i[unlist(cells)], na.rm=T)
  exportSD.S <- sd(rasterLayerHRS.i[unlist(cells)], na.rm=T)
  #floating gC 
  rasterLayerGCF.i <- crop(rasterLayerGC.F,extent(mask))
  cells <- cellFromPolygon(rasterLayerGCF.i, mask[,1])
  gCEstimate.F <- mean(rasterLayerGCF.i[unlist(cells)], na.rm=T)
  gCSD.F <- sd(rasterLayerGCF.i[unlist(cells)], na.rm=T)
  #sinking gC 
  rasterLayerGCS.i <- crop(rasterLayerGC.S,extent(mask))
  cells <- cellFromPolygon(rasterLayerGCS.i, mask[,1])
  gCEstimate.S <- mean(rasterLayerGCS.i[unlist(cells)], na.rm=T)
  gCSD.S <- sd(rasterLayerGCS.i[unlist(cells)], na.rm=T)
  n=as.numeric(lengths(cells))
  #Subtidal.brown.NPP.resampled
  rasterLayerNPP.i <- crop(Subtidal.brown.NPP.resampled,extent(mask))
  cells <- cellFromPolygon(rasterLayerNPP.i, mask[,1])
  NPPEstimate <- mean(rasterLayerNPP.i[unlist(cells)], na.rm=T)
  NPPSD <- sd(rasterLayerNPP.i[unlist(cells)], na.rm=T)
  
  resultsExportEEZOnly <- rbind(resultsExportEEZOnly,data.frame(eez=mask$EEZ[1],country=mask$country[1],
                                                                Basin=mask$seaBasn[1], exportF=exportEstimate.F, 
                                                                exportFSD=exportSD.F, exportS=exportEstimate.S, 
                                                                exportSSD=exportSD.S, gCF=gCEstimate.F, gCF.SD=gCSD.F,
                                                                gCS=gCEstimate.S, gCS.SD=gCSD.S,N=n, NPP=NPPEstimate, NPP.SD=NPPSD))
  
}

#resultsExportEEZOnly

#convert NPP to g
resultsExportEEZOnly$NPP=resultsExportEEZOnly$NPP*1000
resultsExportEEZOnly$NPP.SD=resultsExportEEZOnly$NPP.SD*1000

#calculate the export and gC for combined floating and sinking areas
resultsExportEEZOnly$exportFS=0.5*(resultsExportEEZOnly$exportS+resultsExportEEZOnly$exportF)
resultsExportEEZOnly$gCSF=0.5*(resultsExportEEZOnly$gCS+resultsExportEEZOnly$gCF)

write.csv(resultsExportEEZOnly,file=paste0("resultsExport.EEZ.floatingandsinking.csv"), row.names = FALSE)
save(resultsExportEEZOnly,file=paste0("resultsExport.EEZ.floatingandsinking.RData"))


mean(resultsExportEEZOnly$exportFS,na.rm=T)
mean(resultsExportEEZOnly$gCSF,na.rm=T)


# Ecoregions--------
#do the results by ecoregion

#download from here: https://www.worldwildlife.org/publications/marine-ecoregions-of-the-world-a-bioregionalization-of-coastal-and-shelf-areas
marineEco <- shapefile("Carbon export from seaweed forests to deep marine sinks/Data/MEOW/meow_ecos.shp")
marineEco <- st_as_sf(marineEco)
marineEco <- st_simplify(marineEco, preserveTopology=TRUE, dTolerance = 0.001)
marineEco <- st_make_valid(marineEco)

resultsExportEco2FloatingSinking <- data.frame()

for( i in 1:nrow(marineEco) ) {
  #rasterLayerHR <- rasterLayerHR * raster::area(rasterLayerHR)
  
  cat("\014")
  cat("\n")
  cat("# ---------------------------------\n")
  cat("i: ",i,"out of",nrow(marineEco),"\n")
  cat("# ---------------------------------\n")
  cat("\n")
  
  mask <- marineEco[i,]
  
  #floating export
  rasterLayerHRF.i <- crop(rasterLayerHRF,extent(mask))
  cells <- sfraster::cellFromPolygon(rasterLayerHRF.i, mask[,1])
  exportEstimate.F <- mean(rasterLayerHRF.i[unlist(cells)], na.rm=T)
  exportSD.F <- sd(rasterLayerHRF.i[unlist(cells)], na.rm=T)
  #sinking export
  rasterLayerHRS.i <- crop(rasterLayerHR,extent(mask))
  cells <- sfraster::cellFromPolygon(rasterLayerHRS.i, mask[,1])
  exportEstimate.S <- mean(rasterLayerHRS.i[unlist(cells)], na.rm=T)
  exportSD.S <- sd(rasterLayerHRS.i[unlist(cells)], na.rm=T)
  #floating gC 
  rasterLayerGCF.i <- crop(rasterLayerGC.F,extent(mask))
  cells <- cellFromPolygon(rasterLayerGCF.i, mask[,1])
  gCEstimate.F <- mean(rasterLayerGCF.i[unlist(cells)], na.rm=T)
  gCSD.F <- sd(rasterLayerGCF.i[unlist(cells)], na.rm=T)
  #sinking gC 
  rasterLayerGCS.i <- crop(rasterLayerGC.S,extent(mask))
  cells <- cellFromPolygon(rasterLayerGCS.i, mask[,1])
  gCEstimate.S <- mean(rasterLayerGCS.i[unlist(cells)], na.rm=T)
  gCSD.S <- sd(rasterLayerGCS.i[unlist(cells)], na.rm=T)
  n=as.numeric(lengths(cells))
  #Subtidal.brown.NPP.resampled
  rasterLayerNPP.i <- crop(Subtidal.brown.NPP.resampled,extent(mask))
  cells <- cellFromPolygon(rasterLayerNPP.i, mask[,1])
  NPPEstimate <- mean(rasterLayerNPP.i[unlist(cells)], na.rm=T)
  NPPSD <- sd(rasterLayerNPP.i[unlist(cells)], na.rm=T)
  
  resultsExportEco2FloatingSinking <- rbind(resultsExportEco2FloatingSinking,
    data.frame(marineEco=mask$ECOREGION[1],marineProv=mask$PROVINCE[1],marineRealm=mask$REALM[1],exportF=exportEstimate.F, 
               exportFSD=exportSD.F, exportS=exportEstimate.S, 
               exportSSD=exportSD.S, gCF=gCEstimate.F, gCF.SD=gCSD.F,
               gCS=gCEstimate.S, gCS.SD=gCSD.S,N=n, NPP=NPPEstimate, NPP.SD=NPPSD))
  
}

#convert NPP to g
resultsExportEco2FloatingSinking$NPP=resultsExportEco2FloatingSinking$NPP*1000
resultsExportEco2FloatingSinking$NPP.SD=resultsExportEco2FloatingSinking$NPP.SD*1000

#calculate the export and gC for combined floating and sinking areas
resultsExportEco2FloatingSinking$exportFS=0.5*(resultsExportEco2FloatingSinking$exportS+resultsExportEco2FloatingSinking$exportF)
resultsExportEco2FloatingSinking$gCSF=0.5*(resultsExportEco2FloatingSinking$gCS+resultsExportEco2FloatingSinking$gCF)

#write the results file
write.csv(resultsExportEco2FloatingSinking,file="resultsExportMarEco.floating.sinking.csv", row.names = FALSE)
save(resultsExportEco2FloatingSinking,file=paste0("resultsExportMarEco.floating.sinking.RData"))

mean(resultsExportEco2FloatingSinking$exportFS,na.rm=T)



# -------------------------------------------------#
#### #plot data for specific EEZs  ####               
# -------------------------------------------------#  


rasterLayergc.F <- raster("Outputs/raster.gCexported.2023.10.floating.tif")
rasterLayergc.S <- raster("Outputs/raster.gCexported.2023.10.sinking.tif")

as.data.frame(eez[,2])

#Australia floating
mask <- eez[12,]
plot(mask[2])
rasterLayerHR.i <- crop(rasterLayerF,extent(mask))
plot(rasterLayerHR.i)
rasterLayerHRS.i <- crop(rasterLayer,extent(mask))
plot(rasterLayerHRS.i)
rasterLayerNPP.i <- crop(Subtidal.brown.NPP.resampled,extent(mask))
plot(rasterLayerNPP.i)



#Norway
mask <- marineEco[217,]
plot(mask[2])
rasterLayerHR.i <- crop(rasterLayerHRF,extent(mask))
plot(rasterLayerHR.i)
rasterLayerHRS.i <- crop(rasterLayer,extent(mask))
plot(rasterLayerHRS.i)
rasterLayerNPP.i <- crop(Subtidal.brown.NPP.resampled,extent(mask))
plot(rasterLayerNPP.i)



#Norway EEZ 172
mask <- eez[172,]
plot(mask[2])
rasterLayerHR.i <- crop(rasterLayer,extent(mask))
plot(rasterLayerHR.i)
rasterLayerCRT.i <- crop(rasterLayergc.S,extent(mask))
plot(rasterLayerCRT.i)

# -------------------------------------------------#
#### Create final data set ####               
#### combine export with the seaweed distribution areas from Jorge Assis's model
# -------------------------------------------------#


#load the seaweed areas for EEZs from Jorge for floating, sinking and total area for seaweed species  
resultsAreaCombinedEEZ<-read.csv("Carbon export from seaweed forests to deep marine sinks/resultsAreaCombinedEEZ.csv")
resultsAreaFloatingEEZ<-read.csv("Carbon export from seaweed forests to deep marine sinks/resultsAreaFloatingEEZ.csv")
resultsAreaSinkingEEZ<-read.csv("Carbon export from seaweed forests to deep marine sinks/resultsAreaSinkingEEZ.csv")


#combine all the areas
EEZ.final.ST<-merge(resultsAreaSinkingEEZ, resultsAreaFloatingEEZ, by=c("oceanBasin","eez"))
colnames(EEZ.final.ST)=c('oceanBasin', 'eez', 'area.sinking','area.floating')
EEZ.final<-merge(EEZ.final.ST, resultsAreaCombinedEEZ,by=c("oceanBasin","eez"))


#need an area only estimate, right now floating area includes the floating and sinking combined areas

#calculate the area with both floating and sinking
EEZ.final$floatingandsinking.area<-EEZ.final$area.sinking+EEZ.final$area.floating-EEZ.final$area
EEZ.final$floatingandsinking.area[EEZ.final$floatingandsinking.area<0]=0 #there are some very small negative areas E-14 that should be 0
#sinking only
EEZ.final$sinkingonly.area=EEZ.final$area.sinking-EEZ.final$floatingandsinking.area
EEZ.final$sinkingonly.area[EEZ.final$sinkingonly.area<0]=0
#floating only
EEZ.final$floatingonly.area=EEZ.final$area.floating-EEZ.final$floatingandsinking.area
EEZ.final$floatingonly.area[EEZ.final$floatingonly.area<0]=0

#check totals still match
   #EEZ.final$checktotal=EEZ.final$sinkingonly.area+EEZ.final$floatingonly.area+EEZ.final$floatingandsinking.area
   #cor(EEZ.final$checktotal, EEZ.final$area)
   #they do. cor = 1

#make one sheet just by countries by summing across basins
EEZ.final= EEZ.final[,-1] %>%
  mutate_if(is.numeric, tidyr::replace_na, 0) %>% #in case of having NAs
  group_by(eez) %>%
  summarise_all(., sum, na.rm = TRUE)

EEZ.final=EEZ.final[,c(1,4:7)]

#correct brazil by reducing area to 10%
#the water is so turbid here the seaweed does not extend beyond 3-5 m 
#pers. comm. Alejandro Buschmann, Rita Santos
EEZ.final[EEZ.final$eez=='Brazilian Exclusive Economic Zone',2:5]=EEZ.final[EEZ.final$eez=='Brazilian Exclusive Economic Zone',2:5]*.1

#what is proportion of area with sinking, floating or both species
sinkingonlyproportion=sum(EEZ.final$sinkingonly.area)/sum(EEZ.final$area);sinkingonlyproportion
floatingonlyproportion=sum(EEZ.final$floatingonly.area)/sum(EEZ.final$area);floatingonlyproportion
combinedproportion=sum(EEZ.final$floatingandsinking.area)/sum(EEZ.final$area);combinedproportion
#total area of both species. 46.0%
#total sinking only 19.8%
#total floating only 33.6%

# -------------------------------------------------#
####correct areas for percent rock####               
# -------------------------------------------------#

#load percent rock coastline by nation
percentrockEEZ<-read.csv('resultsExport.EEZ.merged.csv')
percentrockEEZ<-unique(percentrockEEZ[,c(3,5)])

#merge with EEZ export dataset
EEZ.final=merge(EEZ.final, percentrockEEZ, by='eez')

#make the EEZ's without a national percent rock estimate use the average global percent rock
averageprock=mean(EEZ.final$percent.rock,na.rm=T)
EEZ.final$percent.rock[is.na(EEZ.final$percent.rock)]<-averageprock

#correct the different areas for percent rock
ResultsEEZ=EEZ.final
ResultsEEZ$area.corrected.P.rock=ResultsEEZ$area*ResultsEEZ$percent.rock*.01
ResultsEEZ$floatingandsinking.area=EEZ.final$floatingandsinking.area*EEZ.final$percent.rock*.01
ResultsEEZ$sinkingonly.area=EEZ.final$sinkingonly.area*EEZ.final$percent.rock*.01
ResultsEEZ$floatingonly.area=EEZ.final$floatingonly.area*EEZ.final$percent.rock*.01

# -------------------------------------------------#
#####Merge area with export and c G flux datasets #######
# -------------------------------------------------#

#merge area results with export
EEZ<-merge(ResultsEEZ, resultsExportEEZOnly, by=c("eez"))

#calculate single export value for every EEZ using proportional areas for floating and sinking species
EEZ$export<-(EEZ$sinkingonly.area/EEZ$area.corrected.P.rock*EEZ$exportS)+
  (EEZ$floatingonly.area/EEZ$area.corrected.P.rock*EEZ$exportF)+
  (EEZ$floatingandsinking.area/EEZ$area.corrected.P.rock*EEZ$exportFS)

#ensure that the nas did not influence calculation. 
#so places with only floating species will have export from sinking only and vice versa
EEZ$export[EEZ$floatingonly.area==0]=EEZ$exportS[EEZ$floatingonly.area==0]
EEZ$export[EEZ$sinkingonly.area==0]=EEZ$exportF[EEZ$sinkingonly.area==0]

#estimate SD for combined export assumine equal N for both groups
EEZ$exportSF.SD=sqrt((EEZ$exportSSD)^2+(EEZ$exportFSD)^2)


#calculate single SD value for every EEZ 
EEZ$exportcombinedSD<-sqrt(((EEZ$N*EEZ$sinkingonly.area/EEZ$area.corrected.P.rock*(EEZ$exportSSD)^2)+
  ((EEZ$N*EEZ$floatingonly.area/EEZ$area.corrected.P.rock)*EEZ$exportFSD^2)+
  (EEZ$N*EEZ$floatingandsinking.area/EEZ$area.corrected.P.rock*(EEZ$exportSF.SD)^2))/EEZ$N)

EEZ$exportcombinedSD[EEZ$floatingonly.area==0]=EEZ$exportSSD[EEZ$floatingonly.area==0]
EEZ$exportcombinedSD[EEZ$sinkingonly.area==0]=EEZ$exportFSD[EEZ$sinkingonly.area==0]  

EEZ$exportcombinedSD[EEZ$floatingandsinking.area==0]=EEZ$exportSF.SD[EEZ$floatingandsinking.area==0]
#weighted average export

weighted.mean(EEZ$export, EEZ$area.corrected.P.rock, na.rm=T)


#####Summarize total g C export for EEZs ##

#average export for both sinking and floating

EEZ$gCperm2<-(EEZ$sinkingonly.area/EEZ$area.corrected.P.rock*EEZ$gCS)+
  (EEZ$floatingonly.area/EEZ$area.corrected.P.rock*EEZ$gCF)+ 
  (EEZ$floatingandsinking.area/EEZ$area.corrected.P.rock*EEZ$gCSF)

#ensure that the nas did not influence calculation. 
#so places with only floating species will have export from sinking only and vice versa
EEZ$gCperm2[EEZ$floatingonly.area==0]=EEZ$gCS[EEZ$floatingonly.area==0]
EEZ$gCperm2[EEZ$sinkingonly.area==0]=EEZ$gCF[EEZ$sinkingonly.area==0]


#calculate SD for gC estimates
EEZ$gCSF.SD=sqrt((EEZ$gCS.SD)^2+(EEZ$gCF.SD)^2)

#calculate single SD value for every EEZ 
EEZ$gCcombinedSD<-sqrt(((EEZ$N*EEZ$sinkingonly.area/EEZ$area.corrected.P.rock*(EEZ$gCS.SD)^2)+
                          ((EEZ$N*EEZ$floatingonly.area/EEZ$area.corrected.P.rock)*EEZ$gCF.SD^2)+
                          (EEZ$N*EEZ$floatingandsinking.area/EEZ$area.corrected.P.rock*(EEZ$gCSF.SD)^2))/EEZ$N)

#correct points with NAs for some areas
EEZ$gCcombinedSD[EEZ$floatingonly.area==0]=EEZ$gCS.SD[EEZ$floatingonly.area==0]
EEZ$gCcombinedSD[EEZ$sinkingonly.area==0]=EEZ$gCF.SD[EEZ$sinkingonly.area==0]  
EEZ$gCcombinedSD[EEZ$floatingandsinking.area==0]=EEZ$gCSF.SD[EEZ$floatingandsinking.area==0]

weighted.mean(EEZ$gCperm2, EEZ$area.corrected.P.rock, na.rm=T)

# -------------------------------------------------#
#####calculate total TgC by EEZ #######
# -------------------------------------------------#

#convert area km2 to m2
EEZ$gC<-EEZ$sinkingonly.area*EEZ$gCS*1000*1000+(EEZ$floatingonly.area*EEZ$gCF*1000*1000)+
  (EEZ$floatingandsinking.area*EEZ$gCSF*1000*1000)

EEZ$gC[EEZ$eez=='Australian Exclusive Economic Zone']=EEZ$NPP[EEZ$eez=='Australian Exclusive Economic Zone']*EEZ$export[EEZ$eez=='Australian Exclusive Economic Zone']*.01*
  EEZ$area.corrected.P.rock[EEZ$eez=='Australian Exclusive Economic Zone']*1000*1000

#convert g C to Tg C
EEZ$TgC<-EEZ$gC/1E12

sum(EEZ$TgC,na.rm = T)



# -------------------------------------------------#
####Write Results File for EEZ with error###
write.csv(EEZ, 'Outputs/FINALresultsEEZv3.csv')
# -------------------------------------------------#



# -------------------------------------------------#
####  EEZ area and basin  ####               
#this calculates the average export and gC for each EEZ and each ocean basin (n=8). 
#it is broken up into multiple loops so that it will run in quicker stages
#the main reason for dividing by EEZs and basins is larger countries like Australia
#Canada and USA have such large coastlines
#the data are exported into separate results files and then collated
#these calculations are done separately for both floating and sinking species
# -------------------------------------------------#


  
  #floating species basins 1-6
  #southern ocean, south atlantic, south pacific, north pacific, South China and Easter Archipelagic Seas, Indian Ocean
  #saved as resultsExport.EEZ.30m.F.csv
  
  resultsExportF <- data.frame()
  
  for( i in 1:6 ) {
    
    oceanBasins.i <- oceanBasins[i,]
    eez.i <- st_intersection(eez, oceanBasins.i,validate=TRUE)
    
    if( nrow(eez.i) == 0 ) { next }
    
    for( j in 1:nrow(eez.i) ) {
      
      cat("\014")
      cat("\n")
      cat("# ---------------------------------\n")
      cat("i: ",i,"out of",nrow(oceanBasins),"\n")
      cat("j: ",j,"out of",nrow(eez.i),"\n")
      cat("# ---------------------------------\n")
      cat("\n")
      
      mask <- eez.i[j,]
      rasterLayerHR.i <- crop(rasterLayerHRF,extent(mask))
      cells <- cellFromPolygon(rasterLayerHR.i, mask[,1])
      exportEstimate <- mean(rasterLayerHR.i[unlist(cells)], na.rm=T)
      exportSD <- sd(rasterLayerHR.i[unlist(cells)], na.rm=T)
      
      rasterLayerGC.i <- crop(rasterLayerGC.F,extent(mask))
      cells <- cellFromPolygon(rasterLayerGC.i, mask[,1])
      gCEstimate <- mean(rasterLayerGC.i[unlist(cells)], na.rm=T)
      gCSD <- sd(rasterLayerGC.i[unlist(cells)], na.rm=T)
      
      resultsExportF <- rbind(resultsExportF,data.frame(oceanBasin=oceanBasins.i$name[1],eez=mask$EEZ[1],export=exportEstimate, 
                                                        exportSD=exportSD, gC=gCEstimate, gC.SD=gCSD))
      
    }
    
    write.csv(resultsExportF,file=paste0("resultsExport.EEZ.30m.F.csv"), row.names = FALSE)
    save(resultsExportF,file=paste0("resultsExport.EEZ.30m.F.RData"))
  } 
  
  #floating species basins 8-10; Baltic Sea, North Atlantic Ocean, Arctic Ocean
  #saved as resultsExport.EEZ.30m.F2.csv
  
  resultsExportF2 <- data.frame()
  for( i in 8:10 ) {
    
    oceanBasins.i <- oceanBasins[i,]
    eez.i <- st_intersection(eez, oceanBasins.i,validate=TRUE)
    
    if( nrow(eez.i) == 0 ) { next }
    
    for( j in 1:nrow(eez.i) ) {
      
      cat("\014")
      cat("\n")
      cat("# ---------------------------------\n")
      cat("i: ",i,"out of",nrow(oceanBasins),"\n")
      cat("j: ",j,"out of",nrow(eez.i),"\n")
      cat("# ---------------------------------\n")
      cat("\n")
      
      mask <- eez.i[j,]
      rasterLayerHR.i <- crop(rasterLayerHRF,extent(mask))
      cells <- cellFromPolygon(rasterLayerHR.i, mask[,1])
      exportEstimate <- mean(rasterLayerHR.i[unlist(cells)], na.rm=T)
      exportSD <- sd(rasterLayerHR.i[unlist(cells)], na.rm=T)
      
      rasterLayerGC.i <- crop(rasterLayerGC.F,extent(mask))
      cells <- cellFromPolygon(rasterLayerGC.i, mask[,1])
      gCEstimate <- mean(rasterLayerGC.i[unlist(cells)], na.rm=T)
      gCSD <- sd(rasterLayerGC.i[unlist(cells)], na.rm=T)
      
      resultsExportF2 <- rbind(resultsExportF2,data.frame(oceanBasin=oceanBasins.i$name[1],eez=mask$EEZ[1],export=exportEstimate, 
                                                          exportSD=exportSD, gC=gCEstimate, gC.SD=gCSD))
      
    }
    
    write.csv(resultsExportF2,file=paste0("resultsExport.EEZ.30m.F2.v2.csv"), row.names = FALSE)
    save(resultsExportF2,file=paste0("resultsExport.EEZ.30m.F2.v2.RData"))
  }   
  
  #floating species basins 7 Mediterranean region
  #this take a while to run due to the many EEZs and overlapping regions so is done as a separate loop
  #saved as resultsExport.EEZ.30m.F3.csv
  
  resultsExportF3 <- data.frame()
  i=7
  
  oceanBasins.i <- oceanBasins[i,]
  eez.i <- st_intersection(eez, oceanBasins.i,validate=TRUE)
  
  if( nrow(eez.i) == 0 ) { next }
  
  for( j in c(1:2,4:20,28:36)) {
    
    cat("\014")
    cat("\n")
    cat("# ---------------------------------\n")
    cat("i: ",i,"out of",nrow(oceanBasins),"\n")
    cat("j: ",j,"out of",nrow(eez.i),"\n")
    cat("# ---------------------------------\n")
    cat("\n")
    
    mask <- eez.i[j,]
    rasterLayerHR.i <- crop(rasterLayerHRF,extent(mask))
    cells <- cellFromPolygon(rasterLayerHR.i, mask[,1])
    exportEstimate <- mean(rasterLayerHR.i[unlist(cells)], na.rm=T)
    exportSD <- sd(rasterLayerHR.i[unlist(cells)], na.rm=T)
    
    rasterLayerGC.i <- crop(rasterLayerGC.F,extent(mask))
    cells <- cellFromPolygon(rasterLayerGC.i, mask[,1])
    gCEstimate <- mean(rasterLayerGC.i[unlist(cells)], na.rm=T)
    gCSD <- sd(rasterLayerGC.i[unlist(cells)], na.rm=T)
    
    resultsExportF3 <- rbind(resultsExportF3,data.frame(oceanBasin=oceanBasins.i$name[1],eez=mask$EEZ[1],export=exportEstimate, 
                                                        exportSD=exportSD, gC=gCEstimate, gC.SD=gCSD))
    
  }
  
  write.csv(resultsExportF3,file=paste0("resultsExport.EEZ.30m.F3.csv"), row.names = FALSE)
  save(resultsExportF3,file=paste0("resultsExport.EEZ.30m.F3.RData"))
  
  
  
  ##check the difference between average estimates if you use different cut offs of 30 or 50 m - yah! it doesn'tmatter. 
  c= merge(resultsExport, resultsexport50m, by=c('oceanBasin', 'eez'))
  
  plot(c$export.x, c$export.y)
  lm(c$export.x~c$export.y)
  #no difference
  
  plot(c$gC.x, c$gC.y)
  lm(c$gC.x~c$gC.y)
  #no difference
  
  
  #sinking species 
  #basins 1-6 southern ocean, south atlantic, south pacific, north pacific, South China and Easter Archipelagic Seas, Indian Ocean
  #saved as resultsExport.EEZ.30m.csv
  
  resultsExport <- data.frame()
  
  for( i in 1:6 ) {
    
    oceanBasins.i <- oceanBasins[i,]
    eez.i <- st_intersection(eez, oceanBasins.i,validate=TRUE)
    
    if( nrow(eez.i) == 0 ) { next }
    
    for( j in 1:nrow(eez.i) ) {
      
      cat("\014")
      cat("\n")
      cat("# ---------------------------------\n")
      cat("i: ",i,"out of",nrow(oceanBasins),"\n")
      cat("j: ",j,"out of",nrow(eez.i),"\n")
      cat("# ---------------------------------\n")
      cat("\n")
      
      mask <- eez.i[j,]
      rasterLayerHR.i <- crop(rasterLayerHR,extent(mask))
      cells <- cellFromPolygon(rasterLayerHR.i, mask[,1])
      exportEstimate <- mean(rasterLayerHR.i[unlist(cells)], na.rm=T)
      exportSD <- sd(rasterLayerHR.i[unlist(cells)], na.rm=T)
      
      rasterLayerGC.i <- crop(rasterLayerGC.S,extent(mask))
      cells <- cellFromPolygon(rasterLayerGC.i, mask[,1])
      gCEstimate <- mean(rasterLayerGC.i[unlist(cells)], na.rm=T)
      gCSD <- sd(rasterLayerGC.i[unlist(cells)], na.rm=T)
      resultsExport <- rbind(resultsExport,data.frame(oceanBasin=oceanBasins.i$name[1],eez=mask$EEZ[1],export=exportEstimate, 
                                                      exportSD=exportSD, gC=gCEstimate, gC.SD=gCSD))
      
    }
    
    write.csv(resultsExport,file=paste0("resultsExport.EEZ.30m.csv"), row.names = FALSE)
    save(resultsExport,file=paste0("resultsExport.EEZ.30m.RData"))
  }   
  
  
  #sinking species basins 8-10; Baltic Sea, North Atlantic Ocean, Arctic Ocean
  #saved as resultsExport.EEZ2.S.csv
  
  resultsExport2 <- data.frame()
  for( i in 8:10 ) {
    
    oceanBasins.i <- oceanBasins[i,]
    eez.i <- st_intersection(eez, oceanBasins.i,validate=TRUE)
    
    if( nrow(eez.i) == 0 ) { next }
    
    for( j in 1:nrow(eez.i) ) {
      
      cat("\014")
      cat("\n")
      cat("# ---------------------------------\n")
      cat("i: ",i,"out of",nrow(oceanBasins),"\n")
      cat("j: ",j,"out of",nrow(eez.i),"\n")
      cat("# ---------------------------------\n")
      cat("\n")
      
      mask <- eez.i[j,]
      rasterLayerHR.i <- crop(rasterLayerHR,extent(mask))
      cells <- cellFromPolygon(rasterLayerHR.i, mask[,1])
      exportEstimate <- mean(rasterLayerHR.i[unlist(cells)], na.rm=T)
      exportSD <- sd(rasterLayerHR.i[unlist(cells)], na.rm=T)
      rasterLayerGC.i <- crop(rasterLayerGC.S,extent(mask))
      cells <- cellFromPolygon(rasterLayerGC.i, mask[,1])
      gCEstimate <- mean(rasterLayerGC.i[unlist(cells)], na.rm=T)
      gCSD <- sd(rasterLayerGC.i[unlist(cells)], na.rm=T)
      
      resultsExport2 <- rbind(resultsExport2,data.frame(oceanBasin=oceanBasins.i$name[1],eez=mask$EEZ[1],export=exportEstimate, 
                                                        exportSD=exportSD, gC=gCEstimate, gC.SD=gCSD 
                                                        
      ))
      
    }
    
    write.csv(resultsExport2,file=paste0("resultsExport.EEZ2.S.csv"), row.names = FALSE)
    save(resultsExport2,file=paste0("resultsExport.EEZ2.S.RData"))
  }  
  
  #sinking species basins 7 Mediterranean region
  #this take a while to run due to the many EEZs and overlapping regions so is done as a separate loop
  #saved as resultsExport.EEZ3.S.csv
  
  resultsExport3 <- data.frame()
  for( i in 7:7 ) {
    
    oceanBasins.i <- oceanBasins[i,]
    eez.i <- st_intersection(eez, oceanBasins.i,validate=TRUE)
    
    if( nrow(eez.i) == 0 ) { next }
    
    for( j in c(1:2,4:20, 22, 24, 29:36))  {
      
      cat("\014")
      cat("\n")
      cat("# ---------------------------------\n")
      cat("i: ",i,"out of",nrow(oceanBasins),"\n")
      cat("j: ",j,"out of",nrow(eez.i),"\n")
      cat("# ---------------------------------\n")
      cat("\n")
      
      mask <- eez.i[j,]
      rasterLayerHR.i <- crop(rasterLayerHR,extent(mask))
      cells <- cellFromPolygon(rasterLayerHR.i, mask[,1])
      exportEstimate <- mean(rasterLayerHR.i[unlist(cells)], na.rm=T)
      exportSD <- sd(rasterLayerHR.i[unlist(cells)], na.rm=T)
      rasterLayerGC.i <- crop(rasterLayerGC.S,extent(mask))
      cells <- cellFromPolygon(rasterLayerGC.i, mask[,1])
      gCEstimate <- mean(rasterLayerGC.i[unlist(cells)], na.rm=T)
      gCSD <- sd(rasterLayerGC.i[unlist(cells)], na.rm=T)
      #exportquantile <- quantile(rasterLayerHR.i[unlist(cells)],probs=0.9, na.rm=T)
      
      resultsExport3 <- rbind(resultsExport3,data.frame(oceanBasin=oceanBasins.i$name[1],eez=mask$EEZ[1],export=exportEstimate, 
                                                        exportSD=exportSD, gC=gCEstimate, gC.SD=gCSD 
                                                        #    ,area=areaEstimate 
      ))
      
    }
    
    write.csv(resultsExport3,file=paste0("resultsExport.EEZ3.S.csv"), row.names = FALSE)
    save(resultsExport3,file=paste0("resultsExport.EEZ3.S.RData"))
  }   
  
  
  
  
  
##### Merge the datasets for floating and sinking species 
 
  #sinking
  resultsexport50mS<-rbind(resultsExport, resultsExport2,resultsExport3)
  write.csv(resultsexport50mS,file="resultsexport50mS.csv", row.names = FALSE)
  save(resultsexport50mS,file=paste0("resultsexport50mS.RData"))
  
  #floating
  resultsexport50mF<-rbind(resultsExportF, resultsExportF2,resultsExportF3)
  write.csv(resultsexport50mF,file="resultsexport50mF.csv", row.names = FALSE)
  
save(resultsexport50mF,file=paste0("resultsexport50mF.RData"))


# -------------------------------------------------#
#### Create final data set for ecoregions ####               
# -------------------------------------------------#

#### combine export with the seaweed distribution areas from Jorge Assis's model

#summary data by ecoregion: 

##### Estimate the gC by ecoregion
#export file
#resultsExportEco2FloatingSinking
  
#get jorge's area calculations for ecoregions
resultsAreaCombinedMarEco<- read.csv("Carbon export from seaweed forests to deep marine sinks/resultsAreaCombinedMarEco.csv")
resultsAreaSinkingMarEco<- read.csv("Carbon export from seaweed forests to deep marine sinks/resultsAreaSinkingMarEco.csv")
resultsAreaFloatingMarEco<- read.csv("Carbon export from seaweed forests to deep marine sinks/resultsAreaFloatingMarEco.csv")

#merge area files
eco.final.ST<-merge(resultsAreaSinkingMarEco,resultsAreaFloatingMarEco,
                    by=c("marineEco", "marineProv", "marineRealm"))
colnames(eco.final.ST)=c('marineEco', 'marineProv', 'marineRealm','area.sinking', 'area.floating')

eco.final<-merge(eco.final.ST,resultsAreaCombinedMarEco,
                    by=c("marineEco", "marineProv", "marineRealm"))

#need an area only estimate, right now floating area includes the floating and sinking combined areas
#correct for percent rock
#calculate the area with both floating and sinking
eco.final$floatingandsinking.area<-(eco.final$area.sinking+eco.final$area.floating-eco.final$area)*averageprock*.01
#sinking only
eco.final$sinkingonly.area=(eco.final$area.sinking-eco.final$floatingandsinking.area)*averageprock*.01
#floating only
eco.final$floatingonly.area=(eco.final$area.floating-eco.final$floatingandsinking.area)*averageprock*.01

eco.final$area.corrected.P.rock=
  eco.final$floatingandsinking.area+eco.final$sinkingonly.area+eco.final$floatingonly.area

#eco.final$checktotal=eco.final$sinkingonly.area+eco.final$floatingonly.area+eco.final$floatingandsinking.area
#cor(eco.final$checktotal, eco.final$area)
#they do. cor = 1

#fix brazil by adding area correction of 10% reef
eco.final[eco.final$marineProv=='Tropical Southwestern Atlantic',c(4:10)]=
  eco.final[eco.final$marineProv=='Tropical Southwestern Atlantic',c(4:10)]*.1

eco<-merge(eco.final, resultsExportEco2FloatingSinking, by=c("marineEco", "marineProv", "marineRealm"))

eco$export<-(eco$sinkingonly.area/eco$area.corrected.P.rock*eco$exportS)+
  (eco$floatingonly.area/eco$area.corrected.P.rock*eco$exportF)+
  (eco$floatingandsinking.area/eco$area.corrected.P.rock*eco$exportFS)





#ensure that the nas did not influence calculation. 
#so places with only floating species will have export from sinking only and vice versa
eco$export[eco$floatingonly.area==0]=eco$exportS[eco$floatingonly.area==0]
eco$export[eco$sinkingonly.area==0]=eco$exportF[eco$sinkingonly.area==0]


#estimate SD for combined export assumine equal N for both groups
eco$exportSF.SD=sqrt((eco$exportSSD)^2+(eco$exportFSD)^2)

#calculate single SD value for every eco 
eco$exportcombinedSD<-sqrt(((eco$N*eco$sinkingonly.area/eco$area.corrected.P.rock*(eco$exportSSD)^2)+
                              ((eco$N*eco$floatingonly.area/eco$area.corrected.P.rock)*eco$exportFSD^2)+
                              (eco$N*eco$floatingandsinking.area/eco$area.corrected.P.rock*(eco$exportSF.SD)^2))/eco$N)
eco$exportcombinedSD[eco$floatingonly.area==0]=eco$exportSSD[eco$floatingonly.area==0]
eco$exportcombinedSD[eco$sinkingonly.area==0]=eco$exportFSD[eco$sinkingonly.area==0]  


eco$gCperm2<-(eco$sinkingonly.area/eco$area.corrected.P.rock*eco$gCS)+
  (eco$floatingonly.area/eco$area.corrected.P.rock*eco$gCF)+ 
  (eco$floatingandsinking.area/eco$area.corrected.P.rock*eco$gCSF)

eco$gCperm2[eco$floatingonly.area==0]=eco$gCS[eco$floatingonly.area==0]
eco$gCperm2[eco$sinkingonly.area==0]=eco$gCF[eco$sinkingonly.area==0]

# -------------------------------------------------#
#####calculate total TgC by Ecoregion #######
# -------------------------------------------------#

#convert area km2 to m2
eco$gC<-eco$sinkingonly.area*eco$gCS*1000*1000+(eco$floatingonly.area*eco$gCF*1000*1000)+
  (eco$floatingandsinking.area*eco$gCSF*1000*1000)

#convert g C to Tg C
eco$TgC<-eco$gC/1E12
#sum(eco$TgC, na.rm=T)
#cannot sum these because the area corrections are not right because they do not exist for ecoregions. 

# -------------------------------------------------#
####Write Results File for ecoregions ###
write.csv(eco, 'Outputs/FINALresultsEcov3.csv')
# -------------------------------------------------#

#Use the EEZ data for summaries. 
#Summed NPP in TgC for each country
EEZ$NPP.TgC=EEZ$area.corrected.P.rock*1000*1000*EEZ$gCperm2/EEZ$export/0.01/1E12

#global NPP
sum(EEZ$NPP.TgC,na.rm = T)

#global export
sum(EEZ$TgC,na.rm = T)

#average export %
sum(EEZ$TgC,na.rm = T)/sum(EEZ$NPP.TgC,na.rm = T)*100

#the average export (not accounting for spatial differences in NPP) is as follows. 
weighted.mean(EEZ$export, EEZ$area.corrected.P.rock, na.rm=T)

#only the export from currents by floating vs. sinking species
#sinking
weighted.mean(EEZ$exportS, EEZ$area.corrected.P.rock, na.rm=T)
sd(EEZ$exportS, na.rm=T)/sqrt(nrow(EEZ))
#floating
weighted.mean(EEZ$exportF, EEZ$area.corrected.P.rock, na.rm=T)
sd(EEZ$exportF, na.rm=T)/sqrt(nrow(EEZ))

#upper and lower quantiles
quantile(EEZ$export, .1, na.rm=T)
quantile(EEZ$export, .9, na.rm=T)

#the actual gC per m-2

#floating
weighted.mean(EEZ$gCF, (EEZ$floatingandsinking.area*.5)+EEZ$floatingonly.area, na.rm=T)
sd(EEZ$gCF, na.rm=T)/sqrt(nrow(EEZ))

#sinking
weighted.mean(EEZ$gCS, (EEZ$floatingandsinking.area*.5)+EEZ$sinkingonly.area, na.rm=T)
sd(EEZ$gCS, na.rm=T)/sqrt(nrow(EEZ))

#total
weighted.mean(EEZ$gCperm2, EEZ$area.corrected.P.rock, na.rm=T)
sd(EEZ$gCperm2, na.rm=T)/sqrt(nrow(EEZ))

quantile(EEZ$gCperm2, .1, na.rm=T)
quantile(EEZ$gCperm2, .9, na.rm=T)

#what are average export values without accounting for area and NPP
mean(getValues(rasterLayer),na.rm=T)
mean(getValues(rasterLayerF),na.rm=T)
summary(getValues(rasterLayerF))
summary(getValues(rasterLayer))


quantile(getValues(rasterLayer), .10,na.rm=T)
quantile(getValues(rasterLayerF), .10,na.rm=T)


#what is difference between sinking and floating species
#there is no different just based on the export model
weighted.mean(eco$exportF, eco$area.corrected.P.rock, na.rm=T)
sd(eco$exportF,na.rm = T)/sqrt(nrow(eco))#SE

weighted.mean(eco$exportS, eco$area.corrected.P.rock, na.rm=T)
sd(eco$exportS,na.rm = T)/sqrt(nrow(eco))#SE

#but when you weight the average export 
#floating
weighted.mean(eco$exportF, eco$area.corrected.P.rock-eco$sinkingonly.area, na.rm=T)
sd(eco$exportF,na.rm = T)/sqrt(nrow(eco))#SE

weighted.mean(eco$exportS, eco$area.corrected.P.rock-eco$floatingonly.area, na.rm=T)
sd(eco$exportS,na.rm = T)/sqrt(nrow(eco))#SE

#check floating
weighted.mean(eco$exportFS, eco$area.corrected.P.rock-eco$sinkingonly.area, na.rm=T)
sd(eco$exportF,na.rm = T)/sqrt(nrow(eco))#SE

weighted.mean(eco$exportFS, eco$area.corrected.P.rock-eco$floatingonly.area, na.rm=T)
sd(eco$exportS,na.rm = T)/sqrt(nrow(eco))#SE

eco$Proportionfloating = (eco$sinkingonly.area+(eco$floatingandsinking.area*.05))/eco$area.corrected.P.rock


#floating
weighted.mean(EEZ$exportF, EEZ$area.corrected.P.rock-EEZ$sinkingonly.area, na.rm=T)
sd(EEZ$exportF,na.rm = T)/sqrt(nrow(EEZ))#SE

weighted.mean(EEZ$exportS, EEZ$area.corrected.P.rock-EEZ$floatingonly.area, na.rm=T)
sd(EEZ$exportS,na.rm = T)/sqrt(nrow(EEZ))#SE


#what is that in percentages
weighted.mean(EEZ$gCperm2, EEZ$area.corrected.P.rock-EEZ$sinkingonly.area, na.rm=T)/
  weighted.mean(EEZ$NPP, EEZ$area.corrected.P.rock-EEZ$floatingonly.area, na.rm=T)

weighted.mean(EEZ$gCF, EEZ$area.corrected.P.rock-EEZ$sinkingonly.area, na.rm=T)/
  weighted.mean(EEZ$NPP, EEZ$area.corrected.P.rock-EEZ$floatingonly.area, na.rm=T)

weighted.mean(EEZ$gCS, EEZ$area.corrected.P.rock-EEZ$sinkingonly.area, na.rm=T)/
  weighted.mean(EEZ$NPP, EEZ$area.corrected.P.rock-EEZ$floatingonly.area, na.rm=T)



# --------Mesoboundary layer
library(raster)
# load in mesobathypelagic boundary boundary
setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")
meso.data <- read.csv("Datasets/EPI_MESO_BATHY_BOUNDARY.csv")%>% 
  # dplyr::filter(Tax_group=="Brown")%>% #exclude NAs
  glimpse()


meso.fig.export=ggplot() +  
  #geom_polygon(data = meso.data, aes(x=Lon, y=Lat), col='black', size=0.1,fill="grey100")+
  geom_tile(data=meso.data, aes(x=Lon, y=Lat, fill=MESO_BATHY_BOUNDARY_M)) + 
  scale_fill_gradientn(colours = viridis(256, begin= 0.1, end=0.9)) + 
  xlab(expression(paste('longitude (',~degree,')',sep='')))+ylab(expression(paste('latitude (',~degree,')',sep='')))+
  coord_equal() +
  theme()+labs(fill='Sequestration horizon (m)')+
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major = element_blank(),
        plot.background = element_blank(),
        #  legend.position = c(0.5,0.1),legend.direction = 'horizontal',
        legend.position = 'top', legend.key.size = unit(1.2, "cm"), 
        legend.background = element_rect(fill = "gray100", color='black', linetype='solid', size=0.4),
        strip.text = element_text(size = 12),
        legend.key.width=unit(1.2, "cm")) # set length of the unit scale

meso.fig.export


data1 <- data.frame(x = meso.data$Lon,    # Create data frame for raster
                    y = meso.data$Lat,
                    value = meso.data$MESO_BATHY_BOUNDARY_M)
meso.raster <- rasterFromXYZ(data1) 

plot(meso.raster)
#crop to 100 - 400 m depth
values(bathy)[values(bathy < -400)] <- NA
values(bathy)[values(bathy < -50)] <- NA

#
#bathymetryLayerdeep <- raster("Carbon export from seaweed forests to deep marine sinks/GEBCO Bathymetry Global.tif")
#bathymetryLayerdeep[bathymetryLayerdeep < -300] <- NA
#values(bathymetryLayerdeep)[values(bathymetryLayerdeep < -300)] <- NA


rasterLayerSEQH <- resample(meso.raster,bathy, method="ngb")
rasterLayerSEQH <- raster::mask(rasterLayerSEQH,bathy)
plot(rasterLayerSEQH)

resultsSeqH <- data.frame()

#nrow(marineEco)

for( i in 1:3 ) {
  
  #rasterLayerHR <- rasterLayerHR * raster::area(rasterLayerHR)
  
  cat("\014")
  cat("\n")
  cat("# ---------------------------------\n")
  cat("i: ",i,"out of",nrow(marineEco),"\n")
  cat("# ---------------------------------\n")
  cat("\n")

  mask <- marineEco[i,]
  rasterLayerHR.i <- crop(rasterLayerSEQH,extent(mask))
   cells <- sfraster::cellFromPolygon(rasterLayerHR.i, mask[,1])
  meandepth <- mean(rasterLayerHR.i[unlist(cells)], na.rm=T)
  meandepth <- mean(rasterLayerHR.i, na.rm=T)
  depthSD <- sd(rasterLayerHR.i[unlist(cells)], na.rm=T)
  mindepthSD <- min(rasterLayerHR.i[unlist(cells)], na.rm=T)
  
  resultsSeqH <- rbind(resultsSeqH,data.frame(marineEco=mask$ECOREGION[1],marineProv=mask$PROVINCE[1],marineRealm=mask$REALM[1],
                                             SHdepth=meandepth, 
                                             SHdepthSD=depthSD, SHMindepth=mindepthSD))
  
}




# Functions --------
=
 
  #function from gitlab
#I think it is in the fasterize package
  
  cellFromPolygon <- function (object, p, weights = FALSE) {
    which_cells <- function(x) {
      wch <- which(!is.na(raster::values(fasterize::fasterize(x, object))))
      ## this ensures we don't error when no overlap is found
      if (length(wch) < 1L) NULL else wch
    }
    if (inherits(p, "sf")) {
      cells_list <- lapply(split(p, seq_len(nrow(p))), which_cells
      )
      #cellFromPolygon(object, as(p, "Spatial"), weights = weights)
    } else {
      cells_list <- raster::cellFromPolygon(object, p, weights = weights)
    }
    cells_list
  }
  
  