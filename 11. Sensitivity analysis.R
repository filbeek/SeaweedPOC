#sensitivity analysis
setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")

# -------------------------------------------------#
#### run this after running script 06.           
#decomposition sensitivity
# -------------------------------------------------#

crt.raster = brick("Datasets/crt.multiband.tif")
#plot(crt.raster)

# the extent is weird so we need to set it properly 
# this makes the correct projection, same as biooracle rasters
extent(crt.raster) <- extent(-180, 180, -90, 90)

# select only the shallow layers
crt.raster.shallow <- subset(crt.raster, c("ocean_age_200_z_l.2.5" ), drop=F)

# -------------------------------------------------#
#### BOTTOM CURRENTS & BATHYMETRY Sensitivity Analysis               
# -------------------------------------------------#

# Read in bottom current raster from BioOracle
current.raster = raster('Datasets/current.raster.tif')

# Subset the raster by a minimum current threshold
current.subset <- current.raster

#Benthic currents	m/s variance
#Median	0.041
#Average	0.045
#SD	0.012
#SE	0.004
#Min	0.028
#Max	0.067
#N	8

# current.subset[current.raster < 0.062525] <- NA

# Calculate the % areas of seafloor on the shelf without fast enough currents for bedload

bathy.raster=raster('Datasets/bathy.raster.tif')

###### ...Subset by depth ----
## crop by 50 m value. This is matches the seaweed distributions

bathy.subset <- bathy.raster
bathy.subset[bathy.subset<=-50]<-NA ## all the other values give them NA

#this needs a moment to run sometimes... 
crt.resampled= resample(crt.raster.shallow, # the raster you want to fix the resolution of (and keep the values)
                        bathy.subset) #the other raster

crt.resampled # raster is now at the adequate resolution (0.08333)

## now clip the resampled raster for to shallow zone
crt.resampled= mask(crt.resampled, # 
                        bathy.subset) 
#how many areas do we remove by cropping the bottom current out.

#threshold minus SE
x=0.0486-0.004 
current.subset <- current.raster
current.subset[current.raster < x] <- NA
cropped.crt <- mask(crt.resampled, ## large raster first
                    current.subset) # smaller raster second

summary(values(crt.resampled))
summary(values(cropped.crt))
#with lower bound SE NA=9149018

N= 9331200
#((N-3182405)+(N-4611500))/(N-4611500)

((N-9128381)-(N-9149018 ))/(N-9128381)*100
#10.2%

#threshold plus SE
x=0.0486+0.004 
current.subset <- current.raster
current.subset[current.raster < x] <- NA
cropped.crt <- mask(crt.resampled, ## large raster first
                    current.subset) # smaller raster second

summary(values(crt.resampled))
summary(values(cropped.crt))
#with lower bound SE NA=4955923

((N-9128381)-(N-9155638))/(N-9128381)*100
#13.44% zeroed

#average
x=0.0486
current.subset <- current.raster
current.subset[current.raster < x] <- NA
cropped.crt <- mask(crt.resampled, ## large raster first
                    current.subset) # smaller raster second

summary(values(crt.resampled))
summary(values(cropped.crt))
#with lower bound SE NA=9152121
((N-9128381)-(N-9152121))/(N-9128381)*100
#11.7%


####FLOATING LONGEVITY ############################################################
#floating longevity
#sensitivity analysis
#lower SE 9.1
c=31.7-9.1
cropped.crt.depth.floating<-cropped.crt.depth
cropped.crt.depth.floating[cropped.crt.depth.floating <= c] <- c
summary(cropped.crt.depth.floating)
mean(values(cropped.crt.depth), na.rm=T);mean(values(cropped.crt.depth.floating), na.rm=T)


#upper SE 9.1
c=31.7+9.1
cropped.crt.depth.floating<-cropped.crt.depth
cropped.crt.depth.floating[cropped.crt.depth.floating <= c] <- c
summary(cropped.crt.depth.floating)
mean(values(cropped.crt.depth), na.rm=T)-mean(values(cropped.crt.depth.floating), na.rm=T)


#### k values ############################################################
#get x1-x10 from script 06. 
#decomposition rate sensitivity analysis
#only use the upper 5 and lower 5, 5 % quantiles

## we begin by exploring export using single decomposition values, 
## including mean, median and upper and lower 5% quantiles.
raster.export.5q <-((exp(-x2*cropped.crt.depth.floating)))*100*.71
raster.export.median <-((exp(-0.0202*cropped.crt.depth.floating)))*100*.71
raster.export.95q <-((exp(-x9*cropped.crt.depth.floating)))*100*.71
raster.export.mean <-((exp(-0.04548996*cropped.crt.depth.floating)))*100*.71
raster.export.mean <-((exp(-0.04548996*cropped.crt.depth.sinking)))*100*.71

#using mean k + and - SE
XMSEHigher = mean(decomposition.data.brown$k_d, na.rm = T)+(sd(decomposition.data.brown$k_d, na.rm = T)/sqrt(130))
XMSELower = mean(decomposition.data.brown$k_d, na.rm = T)-(sd(decomposition.data.brown$k_d, na.rm = T)/sqrt(130))
raster.export.mean.XMSEHigher <-((exp(-XMSEHigher*cropped.crt.depth.sinking)))*100*.71
raster.export.mean.XMSELower <-((exp(-XMSELower*cropped.crt.depth.sinking)))*100*.71

## compare summary export rates using different decomposition rates
mean(values(raster.export.mean) , na.rm = T)
sd(values(raster.export.mean), na.rm=T)/sqrt(9331200-9.128381e+06) #SE
mean(values(raster.export.mean.XMSELower), na.rm=T)
mean(values(raster.export.mean.XMSEHigher), na.rm=T)
mean(values(raster.export.95q), na.rm=T)
mean(values(raster.export.5q), na.rm=T)

#you should not use these averages to understand average seaweed export because it is not linked to the NPP data or seaweed data yet
originalSinking=C.exported.Raster.S._df 
originalFloating=C.exported.Raster.F._df 

#raster for gC
#raster.export.dist = floating
#raster.export.distnocorrection = sinking

#lower quantiles
raster.export.dist.LOW.F <-((exp(-x1*cropped.crt.depth.floating)*0.2)+(exp(-x2*cropped.crt.depth.floating)*0.2)+(exp(-x3*cropped.crt.depth.floating)*0.2)+
                        (exp(-x4*cropped.crt.depth.floating)*0.2)
                      +(exp(-x5*cropped.crt.depth.floating)*0.2))*100*.71

#upper quantiles
raster.export.dist.HIGH.F <-((exp(-x6*cropped.crt.depth.floating)*0.2)+(exp(-x7*cropped.crt.depth.floating)*0.2)+(exp(-x8*cropped.crt.depth.floating)*0.2)+
                            (exp(-x9*cropped.crt.depth.floating)*0.2)
                          +(exp(-x10*cropped.crt.depth.floating)*0.2))*100*.71


raster.export.dist.LOW.S <-((exp(-x1*cropped.crt.depth.sinking)*0.2)+(exp(-x2*cropped.crt.depth.sinking)*0.2)+(exp(-x3*cropped.crt.depth.sinking)*0.2)+
                              (exp(-x4*cropped.crt.depth.sinking)*0.2)
                            +(exp(-x5*cropped.crt.depth.sinking)*0.2))*100*.71

raster.export.dist.HIGH.S <-((exp(-x6*cropped.crt.depth.sinking)*0.2)+(exp(-x7*cropped.crt.depth.sinking)*0.2)+(exp(-x8*cropped.crt.depth.sinking)*0.2)+
                               (exp(-x9*cropped.crt.depth.sinking)*0.2)
                             +(exp(-x10*cropped.crt.depth.sinking)*0.2))*100*.71

####Compare difference scenarios#
#what is the average difference for sinking
lowerS=mean(values(raster.export.dist.HIGH.S),na.rm=T)
upperS=mean(values(raster.export.dist.LOW.S),na.rm=T)
sinking=mean(values(raster.export.distnocorrection), na.rm=T)

lowerS/sinking
upperS/sinking

lowerF=mean(values(raster.export.dist.HIGH.F),na.rm=T)
upperF=mean(values(raster.export.dist.LOW.F),na.rm=T)
floating=mean(values(raster.export.dist), na.rm=T)

lowerF/floating
upperF/floating


#multiply the export by the NPP, converting the % to a proportion
#convert to g C m-1 so multiply by 1000 then by 0.01

#write the scenario we want to test here: 
#bottom 5 and top 5 q10s for k
testingF=raster.export.dist.LOW.F
testingS=raster.export.dist.LOW.S

#floating
C.exported.Raster.F.test <- overlay(Subtidal.brown.NPP.resampled, testingF,
                               fun=function(r1, r2){return(r1*r2*1000*.01)})

C.exported.Raster.S.test <- overlay(Subtidal.brown.NPP.resampled, testingS,
                               fun=function(r1, r2){return(r1*r2*1000*.01)})


#sinking
C.exported.Raster.S.test_df <- as.data.frame(as(C.exported.Raster.S.test, "SpatialPixelsDataFrame"))
colnames(C.exported.Raster.S.test_df) <- c("C.sequestered", 'lon', 'lat') # name the columns

#get relationship
cor(C.exported.Raster.S._df$C.sequestered,C.exported.Raster.S.test_df$C.sequestered) #highly correlated good! can use the correction 
lm1=lm(C.exported.Raster.S.test_df$C.sequestered~0+C.exported.Raster.S._df$C.sequestered)
lm1
upperexportS=mean(C.exported.Raster.S.test_df$C.sequestered)
exportS= mean(originalSinking$C.sequestered)
upperexportS/exportS

#floating
C.exported.Raster.F.test_df <- as.data.frame(as(C.exported.Raster.F.test, "SpatialPixelsDataFrame"))
colnames(C.exported.Raster.F.test_df) <- c("C.sequestered", 'lon', 'lat') # name the columns

cor(C.exported.Raster.F._df$C.sequestered,C.exported.Raster.F.test_df$C.sequestered)
lm2=lm(C.exported.Raster.F.test_df$C.sequestered~0+C.exported.Raster.F._df$C.sequestered)
lm2

upperexportF=mean(C.exported.Raster.F.test_df$C.sequestered)
exportF=mean(originalFloating$C.sequestered)
upperexportF/exportF

#great mean is very similar to lm 

#now do the lower export scenario of high decomposition
#write the scenario we want to test here: 
testingF=raster.export.dist.HIGH.F
testingS=raster.export.dist.HIGH.S


#floating
C.exported.Raster.F.test <- overlay(Subtidal.brown.NPP.resampled, testingF,
                                    fun=function(r1, r2){return(r1*r2*1000*.01)})

C.exported.Raster.S.test <- overlay(Subtidal.brown.NPP.resampled, testingS,
                                    fun=function(r1, r2){return(r1*r2*1000*.01)})


#sinking
C.exported.Raster.S.test_df <- as.data.frame(as(C.exported.Raster.S.test, "SpatialPixelsDataFrame"))
colnames(C.exported.Raster.S.test_df) <- c("C.sequestered", 'lon', 'lat') # name the columns

lowerexportS=mean(C.exported.Raster.S.test_df$C.sequestered)

#floating
C.exported.Raster.F.test_df <- as.data.frame(as(C.exported.Raster.F.test, "SpatialPixelsDataFrame"))
colnames(C.exported.Raster.F.test_df) <- c("C.sequestered", 'lon', 'lat') # name the columns
lowerexportF= mean(C.exported.Raster.F.test_df$C.sequestered)
lowerexportS/exportS

#upper and lower limit corrections

#get relationship sinking
cor(C.exported.Raster.S._df$C.sequestered,C.exported.Raster.S.test_df$C.sequestered)
lm3=lm(C.exported.Raster.S.test_df$C.sequestered~0+C.exported.Raster.S._df$C.sequestered)
lm3
#summary(lm3)
#this is the worst performing model, but it is still R2 >0.57

cor(C.exported.Raster.F._df$C.sequestered,C.exported.Raster.F.test_df$C.sequestered)
lm4=lm(C.exported.Raster.F.test_df$C.sequestered~0+C.exported.Raster.F._df$C.sequestered)
lm4
#summary(lm4)
lowerexportF=mean(C.exported.Raster.F.test_df$C.sequestered)
lowerexportF/exportF


##############################################################
#Upper and lower 25%

#floating 5, 15, 25%
#lower quantiles 
raster.export.dist.LOW.F <-((exp(-x1*cropped.crt.depth.floating)*1/3)+(exp(-x2*cropped.crt.depth.floating)*1/3)+(exp(-x3*cropped.crt.depth.floating)*1/3))*100*.71

#upper quantiles
raster.export.dist.HIGH.F <-((exp(-x8*cropped.crt.depth.floating)*1/3)+(exp(-x9*cropped.crt.depth.floating)*1/3) +(exp(-x10*cropped.crt.depth.floating)*1/3))*100*.71

#sinking 5, 15, 25%
#lower quantiles 
raster.export.dist.LOW.S <-((exp(-x1*cropped.crt.depth.sinking)*1/3)+(exp(-x2*cropped.crt.depth.sinking)*1/3)+(exp(-x3*cropped.crt.depth.sinking)*1/3))*100*.71
#upper quantiles 
raster.export.dist.HIGH.S <-((exp(-x6*cropped.crt.depth.sinking)*1/3)+(exp(-x7*cropped.crt.depth.sinking)*1/3)+(exp(-x8*cropped.crt.depth.sinking)*1/3))*100*.71

#what is the average difference for sinking
lowerS25=mean(values(raster.export.dist.HIGH.S),na.rm=T)
upperS25=mean(values(raster.export.dist.LOW.S),na.rm=T)
sinking25=mean(values(raster.export.distnocorrection), na.rm=T)

lowerS25/sinking25
upperS25/sinking25

lowerF25=mean(values(raster.export.dist.HIGH.F),na.rm=T)
upperF25=mean(values(raster.export.dist.LOW.F),na.rm=T)
floating25=mean(values(raster.export.dist), na.rm=T)

lowerF25/floating25
upperF25/floating25


#multiply the export by the NPP, converting the % to a proportion
#convert to g C m-1 so multiply by 1000 then by 0.01

#write the scenario we want to test here: 
testingF25=raster.export.dist.LOW.F
testingS25=raster.export.dist.LOW.S


#low k scenario
C.exported.Raster.F.test25 <- overlay(Subtidal.brown.NPP.resampled, testingF25,
                                    fun=function(r1, r2){return(r1*r2*1000*.01)})

C.exported.Raster.S.test25 <- overlay(Subtidal.brown.NPP.resampled, testingS25,
                                    fun=function(r1, r2){return(r1*r2*1000*.01)})

setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")
writeRaster(C.exported.Raster.F.test25, filename = "Outputs/raster.gC.export.floatinglower25.tif", overwrite=TRUE)
writeRaster(C.exported.Raster.S.test25, filename = "Outputs/raster.gC.export.sinkinglower25.tif", overwrite=TRUE)

#sinking
C.exported.Raster.S.test25_df <- as.data.frame(as(C.exported.Raster.S.test25, "SpatialPixelsDataFrame"))
colnames(C.exported.Raster.S.test25_df) <- c("C.sequestered", 'lon', 'lat') # name the columns


#get relationship
cor(C.exported.Raster.S._df$C.sequestered,C.exported.Raster.S.test25_df$C.sequestered)
lm5=lm(C.exported.Raster.S.test25_df$C.sequestered~0+C.exported.Raster.S._df$C.sequestered)
lm5

#floating
C.exported.Raster.F.test25_df <- as.data.frame(as(C.exported.Raster.F.test25, "SpatialPixelsDataFrame"))
colnames(C.exported.Raster.F.test25_df) <- c("C.sequestered", 'lon', 'lat') # name the columns

cor(C.exported.Raster.F._df$C.sequestered,C.exported.Raster.F.test25_df$C.sequestered)
lm6=lm(C.exported.Raster.F.test25_df$C.sequestered~0+C.exported.Raster.F._df$C.sequestered)
lm6

upperexportS25=mean(C.exported.Raster.S.test25_df$C.sequestered)
exportS= mean(originalSinking$C.sequestered)
upperexportS25/exportS


upperexportF25=mean(C.exported.Raster.F.test25_df$C.sequestered)
exportF=mean(originalFloating$C.sequestered)
upperexportF25/exportF

#great mean is very similar to lm 

#now do the lower export scenario of high decomposition
#write the scenario we want to test here: 
#raster.export.dist.HIGH.F
#testingS=raster.export.dist.HIGH.S


#floating
C.exported.Raster.F.test <- overlay(Subtidal.brown.NPP.resampled, raster.export.dist.HIGH.F,
                                    fun=function(r1, r2){return(r1*r2*1000*.01)})

C.exported.Raster.S.test <- overlay(Subtidal.brown.NPP.resampled, raster.export.dist.HIGH.S,
                                    fun=function(r1, r2){return(r1*r2*1000*.01)})

writeRaster(C.exported.Raster.F.test, filename = "Outputs/raster.gC.export.floatingupper25.tif", overwrite=TRUE)
writeRaster(C.exported.Raster.S.test, filename = "Outputs/raster.gC.export.sinkingupper25.tif", overwrite=TRUE)


#get ratios
#sinking
C.exported.Raster.S.test_df <- as.data.frame(as(C.exported.Raster.S.test, "SpatialPixelsDataFrame"))
colnames(C.exported.Raster.S.test_df) <- c("C.sequestered", 'lon', 'lat') # name the columns

lowerexportS25=mean(C.exported.Raster.S.test_df$C.sequestered)

#floating
C.exported.Raster.F.test_df <- as.data.frame(as(C.exported.Raster.F.test, "SpatialPixelsDataFrame"))
colnames(C.exported.Raster.F.test_df) <- c("C.sequestered", 'lon', 'lat') # name the columns

lowerexportF25= mean(C.exported.Raster.F.test_df$C.sequestered)

#upper and lower limit corrections
lowerexportF25/exportF
lowerexportS25/exportS

#get relationship
cor(C.exported.Raster.S._df$C.sequestered,C.exported.Raster.S.test_df$C.sequestered)
lm7=lm(C.exported.Raster.S.test_df$C.sequestered~0+C.exported.Raster.S._df$C.sequestered)
lm7


cor(C.exported.Raster.F._df$C.sequestered,C.exported.Raster.F.test_df$C.sequestered)
lm8=lm(C.exported.Raster.F.test_df$C.sequestered~0+C.exported.Raster.F._df$C.sequestered)
lm8
#summary(lm8)


#stop!
##############################################################
#Upper and lower 15%

#floating 5, 15
#lower quantiles 
raster.export.dist.LOW.F <-((exp(-x1*cropped.crt.depth.floating)*1/3)+(exp(-x2*cropped.crt.depth.floating)*1/3)+(exp(-x3*cropped.crt.depth.floating)*1/3))*100*.71

#upper quantiles
raster.export.dist.HIGH.F <-((exp(-x8*cropped.crt.depth.floating)*1/3)+(exp(-x9*cropped.crt.depth.floating)*1/3) +(exp(-x10*cropped.crt.depth.floating)*1/3))*100*.71

#sinking 5, 15, 25%
#lower quantiles 
raster.export.dist.LOW.S <-((exp(-x1*cropped.crt.depth.sinking)*1/3)+(exp(-x2*cropped.crt.depth.sinking)*1/3)+(exp(-x3*cropped.crt.depth.sinking)*1/3))*100*.71
#upper quantiles 
raster.export.dist.HIGH.S <-((exp(-x6*cropped.crt.depth.sinking)*1/3)+(exp(-x7*cropped.crt.depth.sinking)*1/3)+(exp(-x8*cropped.crt.depth.sinking)*1/3))*100*.71

#what is the average difference for sinking
lowerS25=mean(values(raster.export.dist.HIGH.S),na.rm=T)
upperS25=mean(values(raster.export.dist.LOW.S),na.rm=T)
sinking25=mean(values(raster.export.distnocorrection), na.rm=T)

lowerS25/sinking25
upperS25/sinking25

lowerF25=mean(values(raster.export.dist.HIGH.F),na.rm=T)
upperF25=mean(values(raster.export.dist.LOW.F),na.rm=T)
floating25=mean(values(raster.export.dist), na.rm=T)

lowerF25/floating25
upperF25/floating25


#multiply the export by the NPP, converting the % to a proportion
#convert to g C m-1 so multiply by 1000 then by 0.01

#write the scenario we want to test here: 
testingF25=raster.export.dist.LOW.F
testingS25=raster.export.dist.LOW.S


#low k scenario
C.exported.Raster.F.test25 <- overlay(Subtidal.brown.NPP.resampled, testingF25,
                                      fun=function(r1, r2){return(r1*r2*1000*.01)})

C.exported.Raster.S.test25 <- overlay(Subtidal.brown.NPP.resampled, testingS25,
                                      fun=function(r1, r2){return(r1*r2*1000*.01)})


#sinking
C.exported.Raster.S.test25_df <- as.data.frame(as(C.exported.Raster.S.test25, "SpatialPixelsDataFrame"))
colnames(C.exported.Raster.S.test25_df) <- c("C.sequestered", 'lon', 'lat') # name the columns


#get relationship
cor(C.exported.Raster.S._df$C.sequestered,C.exported.Raster.S.test25_df$C.sequestered)
lm5=lm(C.exported.Raster.S.test25_df$C.sequestered~0+C.exported.Raster.S._df$C.sequestered)
lm5

#floating
C.exported.Raster.F.test25_df <- as.data.frame(as(C.exported.Raster.F.test25, "SpatialPixelsDataFrame"))
colnames(C.exported.Raster.F.test25_df) <- c("C.sequestered", 'lon', 'lat') # name the columns

cor(C.exported.Raster.F._df$C.sequestered,C.exported.Raster.F.test25_df$C.sequestered)
lm6=lm(C.exported.Raster.F.test25_df$C.sequestered~0+C.exported.Raster.F._df$C.sequestered)
lm6

upperexportS25=mean(C.exported.Raster.S.test25_df$C.sequestered)
exportS= mean(originalSinking$C.sequestered)
upperexportS25/exportS


upperexportF25=mean(C.exported.Raster.F.test25_df$C.sequestered)
exportF=mean(originalFloating$C.sequestered)
upperexportF25/exportF

#great mean is very similar to lm 

#now do the lower export scenario of high decomposition
#write the scenario we want to test here: 
testingF=raster.export.dist.HIGH.F
testingS=raster.export.dist.HIGH.S


#floating
C.exported.Raster.F.test <- overlay(Subtidal.brown.NPP.resampled, testingF,
                                    fun=function(r1, r2){return(r1*r2*1000*.01)})

C.exported.Raster.S.test <- overlay(Subtidal.brown.NPP.resampled, testingS,
                                    fun=function(r1, r2){return(r1*r2*1000*.01)})


#sinking
C.exported.Raster.S.test_df <- as.data.frame(as(C.exported.Raster.S.test, "SpatialPixelsDataFrame"))
colnames(C.exported.Raster.S.test_df) <- c("C.sequestered", 'lon', 'lat') # name the columns

lowerexportS25=mean(C.exported.Raster.S.test_df$C.sequestered)

#floating
C.exported.Raster.F.test_df <- as.data.frame(as(C.exported.Raster.F.test, "SpatialPixelsDataFrame"))
colnames(C.exported.Raster.F.test_df) <- c("C.sequestered", 'lon', 'lat') # name the columns

lowerexportF25= mean(C.exported.Raster.F.test_df$C.sequestered)

#upper and lower limit corrections
lowerexportF25/exportF
lowerexportS25/exportS

#get relationship
cor(C.exported.Raster.S._df$C.sequestered,C.exported.Raster.S.test_df$C.sequestered)
lm7=lm(C.exported.Raster.S.test_df$C.sequestered~0+C.exported.Raster.S._df$C.sequestered)
lm7


cor(C.exported.Raster.F._df$C.sequestered,C.exported.Raster.F.test_df$C.sequestered)
lm8=lm(C.exported.Raster.F.test_df$C.sequestered~0+C.exported.Raster.F._df$C.sequestered)
lm8
#summary(lm8)



#Done! 










#GLOBAL CARBON EXPORT ESTIMATES ####               


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

#this completes the decomposition test



###Plot the Results of sensititivity analysis for data
#median decomposition for floating species  

#compare the sinking vs. floating
hist(raster.export.dist.HIGH, main="Size frequency distribution of kelp carbon export", xlim=c(0,100),xlab='% NPP exported across shelf',
     col= "gray50")

hist(raster.export.dist.LOW, main="Size frequency distribution of kelp carbon export", xlim=c(0,100),xlab='% NPP exported across shelf',
     col= "gray50")


#benthic currents
#get df dataframe from 09. 

benthic.test=df

tgC.nocurrents=df$gCF*df$areacorrected.m2/1E12
sum(tgC.nocurrents)

sum(tgC.nocurrents)-sum(df$TgC)


tgC.sinking=df$gCS*df$areacorrected.m2/1E12
sum(tgC.sinking)

sum(tgC.sinking)-sum(df$TgC)
