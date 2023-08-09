######################################################################################
# This script estimates the carbon sink potential of kelp forests 
# It is not yet corrected for spatial area of seaweeds 
# It uses a weighted range of decomposition rates
# It runs without loading 04 and 05 scripts. 
#working directory
setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")

# written by K Filbee-Dexter and A. Pessarrodona  for the EUROMARINE project 21.3.22

#####SET UP###########################################################################

# Load required libraries
library(meowR)
library(tidyverse)
library(sdmpredictors)
library(raster)
library(sf) 
library(ncdf4) # to be able to read netcdfs
library(rgdal)
library(rasterVis)
library(sdmpredictors) #biooracle data layers
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(terra)
library(viridis)


# -------------------------------------------------#
#### NPP & HABITAT ####               
# -------------------------------------------------#

# Read in productivity raster from Pessarrodona et al. 2022
Subtidal.brown.NPP <- raster("Datasets/Subtidal BrownsPresent.tif")
Subtidal.brown.NPP

# Read in brown seaweed forest distribution raster from J. Assis
Subtidal.brown.habitat <- raster(paste0("Datasets/Subtidal Marine forestPresent.tif"))

# -------------------------------------------------#
#### CRT ####               
# -------------------------------------------------#

# Processed coastal residence times from Lui et al. 2019
# Load the raster containing coastal residence times
crt.raster = brick("Datasets/crt.multiband.tif")
# Adjust the extent to match the projection of biooracle rasters
extent(crt.raster) <- extent(-180, 180, -90, 90)
# Select only the shallow layers less then 200 m depth
crt.raster.shallow <- subset(crt.raster, c("ocean_age_200_z_l.2.5"), drop = FALSE)

# -------------------------------------------------#
#### BOTTOM CURRENTS & BATHYMETRY ####               
# -------------------------------------------------#

# Read in bottom current raster from BioOracle
current.raster = raster('Datasets/current.raster.tif')

# Subset the raster by a minimum current threshold
current.subset <- current.raster
current.subset[current.raster < 0.0486] <- NA

# Display summaries of values in the current rasters
summary(getValues(current.subset))


###### ...Subset by depth ----
## crop by 100 m value. This is is coarser then the seaweed distributions 
## but allows proper visualization of the global coastal zone
bathy.raster=raster('Datasets/bathy.raster.tif')
bathy.subset <- bathy.raster
bathy.subset[bathy.subset<=-100]<-NA ## all the other values give them NA
#plot(bathy.subset)


#-----------------------------------------------------------------#
###### CROP BY BATHY & CURRENT RASTER #######
#-----------------------------------------------------------------#
# we now have the CRT raster, but we see that many of the values are for areas offshore where most seaweed wont be found and
#which have lower CRT as they are closer to the shelf. We thus want to crop our CRT raster by a potential seaweed habitat raster,
#we can use the depth as a coarse proxy, and also by areas where bottoms are enough to transport seaweed

#cropped.crt <- mask(crt.raster.shallow, ## large raster first
#                       current.subset) # smaller raster second

# we see that they are not the same resolution, we follow https://gis.stackexchange.com/questions/157557/clip-raster-by-raster-with-data-extraction-and-resolution-change

## first resample to the finer grid
#crt.raster.shallow
#bathy.raster

#this needs a moment to run sometimes... 
crt.resampled= resample(crt.raster.shallow, # the raster you want to fix the resolution of (and keep the values)
                        bathy.raster) #the other raster
#crt.resampled # raster is now at the adequate resolution (0.08333)

#now crop to the depth limit
crt.resampled= mask(crt.resampled, ## large raster first
                        bathy.subset) ## smaller raster second

#plot(crt.resampled)

#look at the areas with low currents
cropped.crt.lowcurrent <- mask(crt.resampled, ## large raster first
                    current.subset) # smaller raster second


#Look at the number of NA cells eliminated
#summary(values(crt.resampled))
#summary(values(cropped.crt.lowcurrent))
(-8973698+9020057)/9020057*100

#-----------------------------------------------------------------#
###### CROP BY A BOX to visualize #######
#-----------------------------------------------------------------#
## lets crop our rasters by the extent of Australia's Great Southern Reef and 
#to have a closer look  

## create a box to crop from
# GSR
min_lon_GSR <- 114
max_lon_GSR <- 155
min_lat_GSR <- -27
max_lat_GSR <- -48

geo_bounds_GSR <- c(left = min_lon_GSR, bottom = min_lat_GSR, right = max_lon_GSR, top = max_lat_GSR)

GSR.grid <- expand.grid(lon_bound_GSR = c(geo_bounds_GSR[1], geo_bounds_GSR[3]), 
                        lat_bound_GSR = c(geo_bounds_GSR[2], geo_bounds_GSR[4]))
coordinates(GSR.grid) <- ~ lon_bound_GSR + lat_bound_GSR # tell R Sites.gridis a spatial object

crt.gsr <- crop(x = crt.resampled, y =  extent(GSR.grid))
plot(crt.gsr)


#-----------------------------------------------------------------#
###### CROP NPP Data #######
#-----------------------------------------------------------------#

#crop the NPP data based on the raw data for AP
setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/NPP, grazing, decomposition")
library(readxl)
productivity<- read_excel("Productivity erosion estimates 11_23_2022.xlsx",1)

productivity=productivity[productivity$Type=='Brown',]
productivity$Avg_ann_prod_kg_C_m2_y<-as.numeric(productivity$Avg_ann_prod_kg_C_m2_y)

#summary data for NPP
#average and lower and upper NPP g C y-1
mean(productivity$Avg_ann_prod_kg_C_m2_y, na.rm=T)*1000
quantile(productivity$Avg_ann_prod_kg_C_m2_y, na.rm=T, .90)*1000
quantile(productivity$Avg_ann_prod_kg_C_m2_y, na.rm=T, .95)*1000

hist(productivity$Avg_ann_prod_kg_C_m2_y*1000)

#95% = 1684 gC m-1 y. We use this to constrain modelled values. 

# Summarize CRT
cropped.crt.depth <- crt.resampled

#plot(cropped.crt.depth)
summary(cropped.crt.depth) 

mean(values(cropped.crt.depth), na.rm=T)
sd(values(cropped.crt.depth), na.rm=T)
N=9331200-8.973698e+06
sd(values(cropped.crt.depth), na.rm=T)/sqrt(N)



#make the cells with NA current data have a large CRT (3000 days) to make export 0

s.new <-
  raster::calc(
    raster::stack(current.subset, cropped.crt.depth),
    fun = function(x)
      if (sum(is.na(x)) > 0)
        x * 100
    else
      x
  )


#now the CRT is large for the slow currents, this compares the two
summary(s.new)

#subset the CRT to just be the CRT raster again
cropped.crt.depth.sinking <- subset(s.new, 2)
summary(cropped.crt.depth.sinking)

#bound the CRT by the maximum CRT (because correction make some values artificially too high)
cropped.crt.depth.sinking[cropped.crt.depth.sinking >= 9.187230e+03] <- 9.187230e+03

#compare, but mean and errors are now artificially high due to the correction
median(values(cropped.crt.depth.sinking), na.rm=T) 
mean(values(cropped.crt.depth.sinking), na.rm=T) #Mean

#average CRT without benthic current correction
median(values(cropped.crt.depth), na.rm=T) 
quantile(values(cropped.crt.depth), na.rm=T, .9)
quantile(values(cropped.crt.depth), na.rm=T, .1)
sd(values(cropped.crt.depth), na.rm=T) #SD


#floating longevity correction
#replace all times shorter than 32 days (floating longevity) with the floating longevity 
cropped.crt.depth.floating<-cropped.crt.depth
cropped.crt.depth.floating[cropped.crt.depth.floating <= 31.7] <- 31.7

#compare, it does not change much
#summary(cropped.crt.depth)
#summary(cropped.crt.depth.floating)

#calculate change in mean CRT from this correction
mean(values(cropped.crt.depth.floating), na.rm=T)-
  mean(values(cropped.crt.depth), na.rm=T)


#now we have two CRT rasters - one for floating and one for sinking species

#look at differences between floating corrected and uncorrected. 
min(values(cropped.crt.depth.floating), na.rm = T);median(values(cropped.crt.depth.floating), na.rm = T)
min(values(cropped.crt.depth.sinking), na.rm = T);median(values(cropped.crt.depth.sinking), na.rm = T)


#crop  NPP by 100 m bathy
Subtidal.brown.NPP.resampled= resample(Subtidal.brown.NPP, # the raster you want to fix the resolution of (and keep the values)
                                       bathy.raster) #the other raster

#Subtidal.brown.NPP.resampled to floating and sinking layers

#remove large NPP values over 1684 gC
#they are from the model and larger than the 95th quantile, l
#likely creating overestimatations in some areas.  

#CRT cropped by the NPP layer

Subtidal.brown.NPP.resampled <- clamp(Subtidal.brown.NPP.resampled, -Inf, 1.684)
max(values(Subtidal.brown.NPP.resampled), na.rm=T)

cropped.crt.NPP <- mask(cropped.crt.depth, ## large raster first
                        Subtidal.brown.NPP.resampled)

cropped.crt.NPP.F <- mask(cropped.crt.depth.floating, ## large raster first
                          Subtidal.brown.NPP.resampled)

cropped.crt.NPP.S <- mask(cropped.crt.depth.sinking, ## large raster first
                          Subtidal.brown.NPP.resampled)

#compare floating versus sinking CRTs. 
median(values(cropped.crt.NPP.S),na.rm=T)
median(values(cropped.crt.NPP.F),na.rm=T)


# -------------------------------------------------#
#### DECOMPOSITION ####               
# -------------------------------------------------#

# Read in decomposition data for brown algae from Pedersen et al.
setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")
decomposition.data <- read.csv("Datasets/Macroalgal decay rates Pedersen Dec 2020_CB-10Jan2022.csv") %>% 
  dplyr::filter(Tax_group == "Brown") %>%  # Exclude NAs
  dplyr::mutate(Genus = stringr::word(Species, 1)) %>% 
  glimpse()

# Summary data for decomposition rates
decomposition.data.summary <- decomposition.data %>% 
  dplyr::filter(!is.na(k_d)) %>%
  group_by(Tax_group) %>%
  dplyr::summarize(
    mean_k_d = mean(k_d),
    median_k_d = median(k_d),
    sd_k_d = sd(k_d),
    observations = length(k_d)
  ) %>% 
  glimpse()

# Capture range of decomposition data using probability distribution
# Dividing the data into 10 categories based on 10 quantiles
decomposition.data.brown = decomposition.data %>% filter(Tax_group == "Brown")
x1 = quantile(decomposition.data.brown$k_d, prob = 0.05, na.rm = TRUE)
x2 = quantile(decomposition.data.brown$k_d, prob = 0.15, na.rm = TRUE)
x3 = quantile(decomposition.data.brown$k_d, prob = 0.25, na.rm = TRUE)
x4 = quantile(decomposition.data.brown$k_d, prob = 0.35, na.rm = TRUE)
x5 = quantile(decomposition.data.brown$k_d, prob = 0.45, na.rm = TRUE)
x6 = quantile(decomposition.data.brown$k_d, prob = 0.55, na.rm = TRUE)
x7 = quantile(decomposition.data.brown$k_d, prob = 0.65, na.rm = TRUE)
x8 = quantile(decomposition.data.brown$k_d, prob = 0.75, na.rm = TRUE)
x9 = quantile(decomposition.data.brown$k_d, prob = 0.85, na.rm = TRUE)
x10 = quantile(decomposition.data.brown$k_d, prob = 0.95, na.rm = TRUE)


#-----------------------------------------------------------------#
###### MULTIPLY RASTER BY DECOMP #######
# this uses decomposition rate and CRT to estimate percent export of NPP
# we do not use an average, but multiply across the weighted range of 10 quantiles decomposition rates for brown algae
# this ensures that we do not obscure the slowly decomposing fraction, which has most change of export. 
# as we are using k values, the equation to get to proportion biomass remaining is exp(-k*days)
# this method underestimates the export because it does not include refractory components
# it is multiplied by 0.71 which is the proportion of NPP released as detritus
# this is the point that the sensitivity analysis is run in 11 script.  
#-----------------------------------------------------------------#



#this calculation accounts for the distribution of k values in 10% quantiles, 
#it is the sum of each weighted proportion of export. 


#every 10 quantile based on CRT for floating seaweeds
raster.export.dist <-((exp(-x1*cropped.crt.depth.floating)*0.1)+(exp(-x2*cropped.crt.depth.floating)*0.1)+(exp(-x3*cropped.crt.depth.floating)*0.1)+(exp(-x4*cropped.crt.depth.floating)*0.1)
                      +(exp(-x5*cropped.crt.depth.floating)*0.1)+(exp(-x6*cropped.crt.depth.floating)*0.1)+(exp(-x7*cropped.crt.depth.floating)*0.1)+(exp(-x8*cropped.crt.depth.floating)*0.1)
                      +(exp(-x9*cropped.crt.depth.floating)*0.1)+(exp(-x10*cropped.crt.depth.floating)*0.1))*100*.71

#every 10 quantile based on CRT for sinking seaweeds
raster.export.distnocorrection <-((exp(-x1*cropped.crt.depth.sinking)*0.1)+(exp(-x2*cropped.crt.depth.sinking)*0.1)+(exp(-x3*cropped.crt.depth.sinking)*0.1)+(exp(-x4*cropped.crt.depth.sinking)*0.1)
                                  +(exp(-x5*cropped.crt.depth.sinking)*0.1)+(exp(-x6*cropped.crt.depth.sinking)*0.1)+(exp(-x7*cropped.crt.depth.sinking)*0.1)+(exp(-x8*cropped.crt.depth.sinking)*0.1)
                                  +(exp(-x9*cropped.crt.depth.sinking)*0.1)+(exp(-x10*cropped.crt.depth.sinking)*0.1))*100*.71


#####Save files
setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")
writeRaster(raster.export.dist, filename = "Outputs/raster.percent.export.dist.floating.2023.11.tif", overwrite=TRUE)
writeRaster(raster.export.distnocorrection, filename = "Outputs/raster.percent.export.dist.sinking.2023.11.tif", overwrite=TRUE)


#compare the export of sinking vs. floating
hist(raster.export.distnocorrection, main="Size frequency distribution of kelp carbon export", xlim=c(0,100),xlab='% NPP exported across shelf',
     col= "gray50")

hist(raster.export.dist, main="Size frequency distribution of kelp carbon export", xlim=c(0,100),xlab='% NPP exported across shelf',
     col= "gray50")


# #summary of global export ----------------------------------------------

# get summary data for export for floating and sinking species

#floating species
raster.export.mean_df <- as.data.frame(as(raster.export.dist, "SpatialPixelsDataFrame"))
colnames(raster.export.mean_df) <- c("export", "lon", "lat") # name the columns
export.summary <- raster.export.mean_df %>% 
  dplyr::summarize(mean_export=mean(export),
                   median_export=median(export),
                   sd_export=sd(export),
                   q10_export=quantile(export,prob=0.1),
                   q90_export=quantile(export,prob=0.9),
                   observations=length(export))%>% 
  glimpse()

raster.export.mean_sinking_df <- as.data.frame(as(raster.export.distnocorrection, "SpatialPixelsDataFrame"))
colnames(raster.export.mean_sinking_df) <- c("export", "lon", "lat") # name the columns
export.summary.sinking <- raster.export.mean_sinking_df %>% 
  dplyr::summarize(mean_export=mean(export),
                   median_export=median(export),
                   sd_export=sd(export),
                   q10_export=quantile(export,prob=0.1),
                   q90_export=quantile(export,prob=0.9),
                   observations=length(export))%>% 
  glimpse()



#-----------------------------------------------------------------#
###### PLOT Global Export for sinking species #######
#-----------------------------------------------------------------#
# load in world shapefile
setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD/Datasets/ne_10m_land")
library(sf)

# Load shapefile
shapename <- read_sf('ne_10m_land.shp')
World.land=shapename

#sinking
#rescale export to better visualize differences
#this now only shows 50% as upper export, instead of 71%
raster.export.dist_rescaled=raster.export.dist
raster.export.dist_rescaled[raster.export.dist_rescaled>=50]<-50 

raster.export.dist_rescaled_df <- as.data.frame(as(raster.export.dist_rescaled, "SpatialPixelsDataFrame"))
colnames(raster.export.dist_rescaled_df) <- c("export", "lon", "lat") # name the columns

#plot %export raw data
World.fig.export=ggplot() +  
  #geom_polygon(data = World.land, aes(x=long, y=lat, group=group), col='black', size=0.1,fill="grey100")+
  geom_tile(data=raster.export.dist_rescaled_df, aes(x=lon, y=lat, fill=export)) + 
  scale_fill_gradientn(colours = viridis(256, begin= 0.1, end=0.9), limits=c(0,50)) + 
  xlab(expression(paste('longitude (',~degree,')',sep='')))+ylab(expression(paste('latitude (',~degree,')',sep='')))+
  coord_equal() +
  theme_classic()+labs(fill='% export')+
  theme_light()+
  theme(legend.position = c(0.1,0.4),legend.direction = 'vertical', legend.key.size = unit(1.2, "cm"), 
        legend.background = element_rect(fill = "gray100", color='black', linetype='solid', size=0.4),
        strip.text = element_text(size = 12),
        legend.key.width=unit(1.2, "cm"))

setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")
tiff('Outputs/World export.floating.tif',res=300, width=24, height=12, units='cm')
World.fig.export
dev.off()



#sinking
raster.export.dist_rescaledS=raster.export.distnocorrection
raster.export.dist_rescaledS[raster.export.dist_rescaledS>=50]<-50 

raster.export.dist_rescaled_df <- as.data.frame(as(raster.export.dist_rescaledS, "SpatialPixelsDataFrame"))
colnames(raster.export.dist_rescaled_df) <- c("export", "lon", "lat") # name the columns

#plot %export raw data
World.fig.exportS=ggplot() +  
  #geom_polygon(data = World.land, aes(x=long, y=lat, group=group), col='black', size=0.1,fill="grey100")+
  geom_tile(data=raster.export.dist_rescaled_df, aes(x=lon, y=lat, fill=export)) + 
  scale_fill_gradientn(colours = viridis(256, begin= 0.1, end=0.9)) + 
  xlab(expression(paste('longitude (',~degree,')',sep='')))+ylab(expression(paste('latitude (',~degree,')',sep='')))+
  coord_equal() +
  labs(fill='% export')+
  theme_light()+
  theme(legend.position = c(0.1,0.4),legend.direction = 'vertical', legend.key.size = unit(1.2, "cm"), 
        legend.background = element_rect(fill = "gray100", color='black', linetype='solid', size=0.4),
        strip.text = element_text(size = 12),
        legend.key.width=unit(1.2, "cm"))

tiff('Outputs/World export.sinking.tif',res=300, width=24, height=12, units='cm')
World.fig.exportS
dev.off()


# -------------------------------------------------#
#### NPP & HABITAT ####               
# -------------------------------------------------#


#productivity raster was read in and resampled to bathy raster above. 

#multiply the export by the NPP, converting the % to a proportion
#convert to g C m-1 so multiply by 1000 then divide by 100

C.exported.Raster.F <- overlay(Subtidal.brown.NPP.resampled, raster.export.distnocorrection,
                               fun=function(r1, r2){return(r1*r2*1000*.01)})

C.exported.Raster.S <- overlay(Subtidal.brown.NPP.resampled, raster.export.dist,
                               fun=function(r1, r2){return(r1*r2*1000*.01)})

plot(C.exported.Raster.S)
plot(C.exported.Raster.F)


#sinking
C.exported.Raster.S._df <- as.data.frame(as(C.exported.Raster.S, "SpatialPixelsDataFrame"))
colnames(C.exported.Raster.S._df) <- c("C.sequestered", 'lon', 'lat') # name the columns

#floating
C.exported.Raster.F._df <- as.data.frame(as(C.exported.Raster.F, "SpatialPixelsDataFrame"))
colnames(C.exported.Raster.F._df) <- c("C.sequestered", 'lon', 'lat') # name the columns

hist(C.exported.Raster.F._df$C.sequestered, col='purple')
mean((C.exported.Raster.F._df$C.sequestered))
sd((C.exported.Raster.F._df$C.sequestered))/sqrt(nrow(C.exported.Raster.F._df))

hist(C.exported.Raster.S._df$C.sequestered, col='purple')
mean((C.exported.Raster.S._df$C.sequestered))
sd((C.exported.Raster.S._df$C.sequestered))/sqrt(nrow(C.exported.Raster.S._df))

# -------------------------------------------------#
#### Save gC files ####               
# -------------------------------------------------#

#output files g C exported per m2 for sinking and floating species
setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")
writeRaster(C.exported.Raster.F, filename = "Outputs/raster.gCexported.2023.10.floating.tif", overwrite=TRUE)
writeRaster(C.exported.Raster.S, filename = "Outputs/raster.gCexported.2023.10.sinking.tif", overwrite=TRUE)


# -------------------------------------------------#
#### GLOBAL CARBON EXPORT ESTIMATES ####               
# -------------------------------------------------#

#what is the carbon store if you use the average exported gC X by the average area in m2 
#using proportion of total area with sinking only, floating only and both types of seaweeds
areasinkingonly=1705227*.12*1000*1000 #12% sinking only
areafloatingonly=1705227*.23*1000*1000 #23% floating only
areaboth=1705227*.65*1000*1000 #65% both

#calculate carbon store in Tg C
carbonstoreSinking=mean((C.exported.Raster.S._df$C.sequestered))*areasinkingonly/1E12 #gC
carbonstoreFloating=mean((C.exported.Raster.F._df$C.sequestered))*areafloatingonly/1E12 #gC
carbonstoreBoth=0.5*(mean((C.exported.Raster.F._df$C.sequestered))+
                       mean((C.exported.Raster.S._df$C.sequestered)))*areaboth/1E12 #gC


carbonstore=carbonstoreSinking+carbonstoreFloating+carbonstoreBoth
carbonstore

