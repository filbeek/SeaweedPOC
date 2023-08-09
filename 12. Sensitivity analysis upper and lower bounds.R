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


#load gC exported layer for floating and sinking species and crop it to bathymetry. 
#upper k and lower export
#rasterLayergc.F <- raster("Outputs/raster.gC.export.floatingupper25.tif")
#rasterLayergc.S <- raster("Outputs/raster.gC.export.sinkingupper25.tif")

#lower k and upper export
rasterLayergc.F <- raster("Outputs/raster.gC.export.floatinglower25.tif")
rasterLayergc.S <- raster("Outputs/raster.gC.export.sinkinglower25.tif")

rasterLayerGC.S <- resample(rasterLayergc.S,bathymetryLayer, method="ngb")
rasterLayerGC.S <- raster::mask(rasterLayerGC.S,bathymetryLayer)

rasterLayerGC.F <- resample(rasterLayergc.F,bathymetryLayer, method="ngb")
rasterLayerGC.F <- raster::mask(rasterLayerGC.F,bathymetryLayer)

mean(values(rasterLayerGC.F), na.rm = T)
mean(values(rasterLayerGC.S), na.rm = T)

# load shapefiles to get boundaries for ocean basins and the worlds EEZs--------
# no need to run every time
# this can take a while to run
eez <- shapefile("EEZ/EEZGlobal.shp")

sf::sf_use_s2(FALSE)


eez <- st_as_sf(eez)
eez <- st_simplify(eez, preserveTopology=TRUE, dTolerance = 0.001)
eez <- st_make_valid(eez)

# -------------------------------------------------#
#### Calculate export for each EEZ area  ####               
#this calculates the average percent export and gC exported for each EEZ WITHOUT dividing by ocean basin 
#it only runs EEZs with seaweed areas. 1:23,25:174,178:186, 191:278
#these calculations are done separately for both floating and sinking species
# -------------------------------------------------#

#floating species 


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
  
  resultsExportEEZOnly <- rbind(resultsExportEEZOnly,data.frame(eez=mask$EEZ[1],country=mask$country[1],
                                                                Basin=mask$seaBasn[1], gCF=gCEstimate.F, gCF.SD=gCSD.F,
                                                                gCS=gCEstimate.S, gCS.SD=gCSD.S,N=n))
  
}

#resultsExportEEZOnly

#convert NPP to g

#calculate the export and gC for combined floating and sinking areas
resultsExportEEZOnly$gCSF=0.5*(resultsExportEEZOnly$gCS+resultsExportEEZOnly$gCF)

write.csv(resultsExportEEZOnly,file=paste0("resultsExport.EEZ.floatingandsinkingUpper25.csv"), row.names = FALSE)


#upper k and lower export
rasterLayergc.F <- raster("Outputs/raster.gC.export.floatingupper25.tif")
rasterLayergc.S <- raster("Outputs/raster.gC.export.sinkingupper25.tif")


rasterLayerGC.S <- resample(rasterLayergc.S,bathymetryLayer, method="ngb")
rasterLayerGC.S <- raster::mask(rasterLayerGC.S,bathymetryLayer)

rasterLayerGC.F <- resample(rasterLayergc.F,bathymetryLayer, method="ngb")
rasterLayerGC.F <- raster::mask(rasterLayerGC.F,bathymetryLayer)

mean(values(rasterLayerGC.F), na.rm = T)
mean(values(rasterLayerGC.S), na.rm = T)



# -------------------------------------------------#
#### Calculate export for each EEZ area  ####               
#this calculates the average percent export and gC exported for each EEZ WITHOUT dividing by ocean basin 
#it only runs EEZs with seaweed areas. 1:23,25:174,178:186, 191:278
#these calculations are done separately for both floating and sinking species
# -------------------------------------------------#

#floating species 


resultsExportEEZOnlyLower= data.frame()

for( i in c(1:23,25:174,178:186, 191:278) ) {
  
  cat("\014")
  cat("\n")
  cat("# ---------------------------------\n")
  cat("i: ",i,"out of",nrow(eez),"\n")
  cat("# ---------------------------------\n")
  cat("\n")
  
  #i=32
  mask <- eez[i,]
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
  
  resultsExportEEZOnlyLower <- rbind(resultsExportEEZOnlyLower,data.frame(eez=mask$EEZ[1],country=mask$country[1],
                                                                Basin=mask$seaBasn[1], gCF=gCEstimate.F, gCF.SD=gCSD.F,
                                                                gCS=gCEstimate.S, gCS.SD=gCSD.S,N=n))
  
}

#resultsExportEEZOnlyLower



#calculate the export and gC for combined floating and sinking areas
resultsExportEEZOnlyLower$gCSF=0.5*(resultsExportEEZOnlyLower$gCS+resultsExportEEZOnlyLower$gCF)

x=merge(resultsExportEEZOnlyLower, resultsExportEEZOnly, by=c('eez', 'country', 'Basin', 'N'))

mean(x$gCF.x, na.rm=T)
mean(x$gCF.y, na.rm=T)

write.csv(x,file=paste0("resultsExport.EEZ.floatingandsinkingLowerUpper25.csv"), row.names = FALSE)


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

#sinking only
EEZ.final$sinkingonly.area=EEZ.final$area.sinking-EEZ.final$floatingandsinking.area

#floating only
EEZ.final$floatingonly.area=EEZ.final$area.floating-EEZ.final$floatingandsinking.area

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
EEZ<-merge(ResultsEEZ, x, by=c("eez"))

#####Summarize total g C export for EEZs ##

#average export for both sinking and floating

EEZ$gCperm2lower25<-(EEZ$sinkingonly.area/EEZ$area.corrected.P.rock*EEZ$gCS.x)+
  (EEZ$floatingonly.area/EEZ$area.corrected.P.rock*EEZ$gCF.x)+ 
  (EEZ$floatingandsinking.area/EEZ$area.corrected.P.rock*EEZ$gCSF.x)

#ensure that the nas did not influence calculation. 
#so places with only floating species will have export from sinking only and vice versa
EEZ$gCperm2lower25[EEZ$floatingonly.area==0]=EEZ$gCS.x[EEZ$floatingonly.area==0]
EEZ$gCperm2gCperm2lower25[EEZ$sinkingonly.area==0]=EEZ$gCF.x[EEZ$sinkingonly.area==0]


#calculate SD for gC estimates
EEZ$gCSF.SD.x=sqrt((EEZ$gCS.SD.x)^2+(EEZ$gCF.SD.x)^2)

#calculate single SD value for every EEZ 
EEZ$gCcombinedSD.x<-sqrt(((EEZ$N*EEZ$sinkingonly.area/EEZ$area.corrected.P.rock*(EEZ$gCS.SD.x)^2)+
                          ((EEZ$N*EEZ$floatingonly.area/EEZ$area.corrected.P.rock)*EEZ$gCF.SD.x^2)+
                          (EEZ$N*EEZ$floatingandsinking.area/EEZ$area.corrected.P.rock*(EEZ$gCSF.SD.x)^2))/EEZ$N)

#correct points with NAs for some areas
EEZ$gCcombinedSD.x[EEZ$floatingonly.area==0]=EEZ$gCS.SD.x[EEZ$floatingonly.area==0]
EEZ$gCcombinedSD.x[EEZ$sinkingonly.area==0]=EEZ$gCF.SD.x[EEZ$sinkingonly.area==0]  
EEZ$gCcombinedSD.x[EEZ$floatingandsinking.area==0]=EEZ$gCSF.SD.x[EEZ$floatingandsinking.area==0]





#upper estimate average export for both sinking and floating

EEZ$gCperm2upper25<-(EEZ$sinkingonly.area/EEZ$area.corrected.P.rock*EEZ$gCS.y)+
  (EEZ$floatingonly.area/EEZ$area.corrected.P.rock*EEZ$gCF.y)+ 
  (EEZ$floatingandsinking.area/EEZ$area.corrected.P.rock*EEZ$gCSF.y)

#ensure that the nas did not influence calculation. 
#so places with only floating species will have export from sinking only and vice versa
EEZ$gCperm2upper25[EEZ$floatingonly.area==0]=EEZ$gCS.y[EEZ$floatingonly.area==0]
EEZ$gCperm2upper25[EEZ$sinkingonly.area==0]=EEZ$gCF.y[EEZ$sinkingonly.area==0]


#calculate SD for gC estimates
EEZ$gCSF.SD.y=sqrt((EEZ$gCS.SD.y)^2+(EEZ$gCF.SD.y)^2)

#calculate single SD value for every EEZ 
EEZ$gCcombinedSD.y<-sqrt(((EEZ$N*EEZ$sinkingonly.area/EEZ$area.corrected.P.rock*(EEZ$gCS.SD.y)^2)+
                            ((EEZ$N*EEZ$floatingonly.area/EEZ$area.corrected.P.rock)*EEZ$gCF.SD.y^2)+
                            (EEZ$N*EEZ$floatingandsinking.area/EEZ$area.corrected.P.rock*(EEZ$gCSF.SD.y)^2))/EEZ$N)

#correct points with NAs for some areas
EEZ$gCcombinedSD.y[EEZ$floatingonly.area==0]=EEZ$gCS.SD.y[EEZ$floatingonly.area==0]
EEZ$gCcombinedSD.y[EEZ$sinkingonly.area==0]=EEZ$gCF.SD.y[EEZ$sinkingonly.area==0]  
EEZ$gCcombinedSD.y[EEZ$floatingandsinking.area==0]=EEZ$gCSF.SD.y[EEZ$floatingandsinking.area==0]


# -------------------------------------------------#
#####calculate total TgC by EEZ #######
# -------------------------------------------------#


#convert area km2 to m2
EEZ$gC.lower25<-EEZ$sinkingonly.area*EEZ$gCS.x*1000*1000+(EEZ$floatingonly.area*EEZ$gCF.x*1000*1000)+
  (EEZ$floatingandsinking.area*EEZ$gCSF.x*1000*1000)
#convert area km2 to m2
EEZ$gC.upper25<-EEZ$sinkingonly.area*EEZ$gCS.y*1000*1000+(EEZ$floatingonly.area*EEZ$gCF.y*1000*1000)+
  (EEZ$floatingandsinking.area*EEZ$gCSF.y*1000*1000)

#convert g C to Tg C
EEZ$TgC.lower25<-EEZ$gC.lower25/1E12
sum(EEZ$TgC.lower25,na.rm = T)


#convert g C to Tg C
EEZ$TgC.upper25<-EEZ$gC.upper25/1E12

sum(EEZ$TgC.upper25,na.rm = T)




# -------------------------------------------------#
####Write Results File for EEZ with error###
write.csv(EEZ, 'Outputs/FINALresultsEEZv3lowerandupper.csv')
# -------------------------------------------------#

#STOP


# -------------------------------------------------#
#### Create final data set for ecoregions ####               
# -------------------------------------------------#

#### combine export with the seaweed distribution areas from Jorge Assis's model

#summary data by ecoregion: 

##### Estimate the gC by ecoregion
#export file
resultsExportEco2FloatingSinking

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

