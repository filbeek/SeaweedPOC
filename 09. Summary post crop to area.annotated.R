######################################################################################
# this script maps the results for export by nation and ecoregion
# it includes estimates of ventilation rates at 208 and 530 m depth
# written by K Filbee-Dexter and A. Pessarrodona  for the EUROMARINE project 21.3.22

######################################################################################


setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")
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
library(sf)
library(colorspace)
library(spatstat.utils)
library(ggmap)
library(metR) 

library(maps)
library(ggplot2)
library(sf)        # for manipulation of simple features objects
library(usmap)
library(rnaturalearth)

#load EEZ results
df=read.csv ('Outputs/FINALresultsEEZv3.csv')
sum(df$area.corrected.P.rock,na.rm = T)
#1.9 Million km2

#uncertainity
#add upper and lower errors using decomposition rate uncertainty estimates from sensitivity analysis. 
upperlower=read.csv('Outputs/FINALresultsEEZv3lowerandupper.csv')

upperlower=upperlower[,c('eez','TgC.lower25','TgC.upper25')]

df=merge(df,upperlower, by='eez')


#difference in total
sum(df$TgC.upper25, na.rm=T)
sum(df$TgC.lower25, na.rm=T)

#now add the SE for % NPP released as detritus
df$TgC.upperwNPP=df$TgC.upper25*(97.6)/71
df$TgC.lowerwNPP=df$TgC.lower25*(50)/71

sum(df$TgC.upperwNPP, na.rm=T)
sum(df$TgC.lowerwNPP, na.rm=T)

#get weighted means
df=df[!is.na(df$TgC),]
#total export
weighted.mean(df$export, df$area.corrected.P.rock, na.rm=T)
sd(df$export, na.rm=T)/sqrt(nrow(df))
quantile(df$export, .9);quantile(df$export, .1)

#compare floating to sinking
weighted.mean(df$exportF, df$area.corrected.P.rock, na.rm=T)
sd(df$exportF, na.rm=T)/sqrt(nrow(df))
weighted.mean(df$exportS, df$area.corrected.P.rock, na.rm=T)
sd(df$exportS, na.rm=T)/sqrt(nrow(df))

weighted.mean(df$gCperm2, df$area.corrected.P.rock, na.rm=T)
sd(df$gCperm2, na.rm=T)/sqrt(nrow(df))
quantile(df$gCperm2, .9);quantile(df$gCperm2, .1)

#compare floating to sinking
weighted.mean(df$gCF, df$floatingandsinking.area*.5+df$floatingonly.area, na.rm=T)
sd(df$gCF, na.rm=T)/sqrt(nrow(df))
weighted.mean(df$gCS, df$floatingandsinking.area*.5+df$sinkingonly.area, na.rm=T)
sd(df$gCS, na.rm=T)/sqrt(nrow(df))

write.csv(df, 'Outputs/Exportbynationwithuncertainty.csv')

NPPglobal=sum(df$NPP*df$area.corrected.P.rock*1000*1000/1E12)
globalCsq/NPPglobal

df$areacorrected.m2<-df$area.corrected.P.rock*1000*1000

globalCsq=sum(df$TgC, na.rm = T) #Tg

#C seq over the last 20 years
globalCsq*20

#in tonnes
globalCsq*1000000
#in gigatonns
globalCsq*.001

#percentage of the biological pump
globalCsq*1000000/11E9*100  #11 Gt C yr−1 billion tons Basa 2018 https://doi.org/10.3390/su10030869
globalCsq*.001/11*100

#updated BC numbers
bluecarbon=41+12.63+35.31
globalCsq+bluecarbon


#load sequestration horizon results
seq<-read_excel("Outputs/sequestration.horizon.ecoregion.xlsx" )
eco2=read.csv ('Outputs/FINALresultsEcov3.csv')

max(eco2$NPP,na.rm =T)


# Combine seq estimates with export --------------------------------------

#merge export by ecoregion with the area calculation (the area columns in the resultsexport file are wrong) 
eco2=merge(eco2, seq,  by=c('marineEco', 'marineRealm', 'marineProv'))



#calculate the percent C reaching different seq places

eco2$seq.208m=as.factor(eco2$seq.208m.y)
eco2$seq.530m=as.factor(eco2$seq.530m.y)

seq208=eco2 %>% 
  dplyr::group_by(seq.208m) %>%
  dplyr:: summarise(Tgseq208m=sum(TgC, na.rm = T), N208=n())


seq530=eco2 %>% 
  dplyr::group_by(seq.530m) %>%
  dplyr:: summarise(Tgseq530m=sum(TgC, na.rm = T), N530=n())

#make a merged dataset
eco3=eco2
eco2$Seqdepth=208
eco2$seq.years=eco2$seq.208m.y

eco3$Seqdepth=530
eco3$seq.years=eco3$seq.530m.y
mergedseq=rbind(eco2, eco3)

mergedseq$Seqdepth=as.factor(mergedseq$Seqdepth)

#fix white bars in plot by adding in dummy rows with zero seaweed carbon values for certain years
x=nrow(mergedseq)
mergedseq[x+1,]=mergedseq[x,]
mergedseq$Seqdepth[x+1]='208'
mergedseq$exportgC[x+1]=0
mergedseq$seq.years[x+1]=181.25

mergedseq[x+2,]=mergedseq[x,]
mergedseq$Seqdepth[x+2]='208'
mergedseq$exportgC[x+2]=0
mergedseq$seq.years[x+2]=218.75

mergedseq[x+3,]=mergedseq[x,]
mergedseq$Seqdepth[x+3]='530'
mergedseq$exportgC[x+3]=0
mergedseq$seq.years[x+3]=243.75

#eco$exportgC<-eco$gCperm2
correction=sum(df$TgC,na.rm=T)/sum(eco2$TgC,na.rm=T)

#Percent carbon reaching areas beyond sequestration horizon
#total over 25 years
sum(eco2$TgC[eco2$seq.208m.y>=25], na.rm=T)/sum(eco2$TgC, na.rm=T)*100
#over 100 years
sum(eco2$TgC[eco2$seq.208m.y>=100], na.rm=T)/sum(eco2$TgC, na.rm=T)*100


#over 25 years for 530
sum(eco2$TgC[eco2$seq.530m.y>=25], na.rm=T)/sum(eco2$TgC, na.rm=T)*100
#over 100 years
sum(eco2$TgC[eco2$seq.530m.y>=100], na.rm=T)/sum(eco2$TgC, na.rm=T)*100


#group by marine realm to see range of seq times
seq.summary=eco2 %>% 
  dplyr::group_by(marineRealm) %>%
  dplyr:: summarise(average.seq208=mean(seq.208m.y, na.rm = T),
                    SD.seq208=sd(seq.208m.y, na.rm = T),
                    average.seq530=mean(seq.530m.y, na.rm = T),
                    SD.seq530=sd(seq.530m.y, na.rm = T),N=n())

seq.summary$SE.208=seq.summary$SD.seq208/sqrt(seq.summary$N)
seq.summary$SE.530=seq.summary$SD.seq530/sqrt(seq.summary$N)
seq.summary

#make cumulative curve
seq.summary=mergedseq %>% 
  dplyr::group_by(seq.years, Seqdepth) %>%
  dplyr:: summarise(gCseq=sum(TgC, na.rm = T)*correction)

seq.summary=seq.summary[!is.na(seq.summary$seq.years),]

#calculate the TgC lasting over 100 years at both depths. between 93 and 106 m
s530=sum(seq.summary[seq.summary$seq.years>=100&seq.summary$Seqdepth=='530',3])+
  sum(seq.summary[seq.summary$seq.years==93.75&seq.summary$Seqdepth=='530',3])*.5
s530
s530/56*100

s208=(sum(seq.summary[seq.summary$seq.years==93.75&seq.summary$Seqdepth=='208',3])*.5+
sum(seq.summary[seq.summary$seq.years>=100&seq.summary$Seqdepth=='208',3]))
s208

s208/56*100
#between 8 and 31.

s530=(sum(seq.summary[seq.summary$seq.years==18.75&seq.summary$Seqdepth=='530',3])*.5+
        sum(seq.summary[seq.summary$seq.years>=31.25&seq.summary$Seqdepth=='530',3]))
s530
s530/56*100

s208=(sum(seq.summary[seq.summary$seq.years==18.75&seq.summary$Seqdepth=='208',3])*.5+
        sum(seq.summary[seq.summary$seq.years>=31.25&seq.summary$Seqdepth=='208',3]))
s208

s208/56*100



# -------------------------------------------------#
#### plot the cumulative sequestration ####               
# -------------------------------------------------#


cumSeq=ggplot() + 
  geom_line(col='#44AA99', linewidth=1.5, data=seq.summary[seq.summary$Seqdepth=='208',], 
            aes(x=seq.years, y= revcumsum(gCseq))) +
  geom_line(col='#332288', linewidth=1.5,data=seq.summary[seq.summary$Seqdepth=='530',],
            aes(x=seq.years, y= revcumsum(gCseq))) +
  ylab(expression(Tg~seaweed~C~y^-1~exported))+
  xlab ('median carbon sequestration timescales (years)')+
  geom_vline(aes(xintercept =100), col='orangered', size=1)+
  geom_vline(aes(xintercept =25), col='orangered')+
  scale_x_continuous(limits = c(0,255), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()

cumSeq


#order
mergedseq$marineRealm <- factor(mergedseq$marineRealm, 
                                  levels = c("Temperate Southern Africa", "Eastern Indo-Pacific",   "Tropical Atlantic", 
                                             "Central Indo-Pacific" ,  "Tropical Eastern Pacific" ,  "Temperate Australasia",  
                                             "Temperate Northern Pacific", "Temperate Northern Atlantic",   "Temperate South America" ,    
                                             "Western Indo-Pacific","Southern Ocean" ,   "Arctic"  ))
unique(mergedseq$marineRealm)

realmseq=ggplot(mergedseq, aes(x=marineRealm, y= seq.years, fill=Seqdepth)) + 
  geom_boxplot() +theme_classic() +ylab ('median carbon sequestration timescales (years)')+xlab("")+
  scale_fill_manual(values=c( '#44AA99','#332288'))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,255)) + coord_flip()+
  theme_classic()+labs(fill="Depth of remineralization (m): ")

cumSeq.fig= ggarrange( cumSeq, realmseq,nrow=1, labels = c('A', 'B'), common.legend = T, widths = c(0.8, 1.2))

tiff('Outputs/cumulative.seq.time.tif', width = 24, height = 12, units = 'cm', res=300)
cumSeq.fig
dev.off()

###########################################################################################
#read summary table by EEZ and percent rock

TotalNPP=df$NPP*df$area.corrected.P.rock*1000*1000
sum(TotalNPP,na.rm=T)/1E12

# Refactory POC--------------------------------------
#estimate refactory of 5% TgC
sum(TotalNPP)*0.05/1E12

#estimate refactory of 1% TgC
sum(TotalNPP)*0.01/1E12

#df$area.corrected.P.rock.test=df$floatingandsinking.area+df$floatingonly.area+df$sinkingonly.area
#cor(df$area.corrected.P.rock.test,df$area.corrected.P.rock)


# -------------------------------------------------#
#### map EXPORT by ecoregion ####               
# -------------------------------------------------#


setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD/Datasets/ne_10m_land")

shapename <- read_sf('ne_10m_land.shp')
World.land=shapename

#% exported
#combine data with ecoregions polygon file
setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD/Carbon export from seaweed forests to deep marine sinks/Data/MEOW")
regions <- readOGR("meow_ecos.shp")
regions@data$id = rownames(regions@data)
regions.points = fortify(regions, ECOREGION="id")
regions.df = join(regions.points, regions@data, by="id")
regions.df3 <- merge(regions.df, eco2, by.x='ECOREGION', by.y='marineEco', .drop=T)
regions.df2=regions.df3#make a copy

#plot export
cols=rev(brewer.pal(11, "Spectral")[3:10])
min=min(regions.df3$export, na.rm=T)

#make maximums similar color for better visualization
regions.df3= regions.df3[!is.na(regions.df3$export),]
regions.df3$export[regions.df3$export>=50]<-51

map.export=ggplot() +  
  geom_polygon(data=regions.df3%>%filter(regions.df3$NPP>0), mapping=aes(long,lat,group=group,fill=export)) +
  geom_path(regions.df3, mapping=aes(long,lat,group=group), col='gray100')+
  scale_fill_gradientn(colors = cols, name = "% NPP exported >200 m", breaks=c(min, 10,20,30,40,50), labels=c('0','10', '20', '30', '40','>50' ))+
  scale_x_longitude(breaks = seq(-180,180,45))+
  scale_y_latitude(breaks = seq(-80,80,40))+
  xlab(expression(paste('longitude (',~degree,')',sep='')))+ylab(expression(paste('latitude (',~degree,')',sep='')))+
  #ylim(-90,88)+
  geom_sf(data=World.land, col='black', size=0.1,fill="grey100")+
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major = element_blank(),
        plot.background = element_blank(), legend.position = c(0.5,1.1),legend.direction = "horizontal", plot.margin = unit(c(2,.1,.1,.1), "cm"), 
        legend.margin=margin(t = 0, unit='cm'), legend.title = element_text(size = 12),legend.text = element_text(size = 12), 
        legend.key.width =  unit(2,'cm'))

map.export

setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")

tiff('Outputs/map.export.percent.tif', width = 30, height = 15, units = 'cm', res=300)
map.export
dev.off()



#crop exported carbon to 90th percentile for better visualization
regions.df3=regions.df2
regions.df3$gCperm2[regions.df3$gCperm2>=501]=501

df3= regions.df3[!is.na(regions.df3$gCperm2),]

#total carbon
map.gC=ggplot() +  
  geom_polygon(data=df3%>%filter(df3$NPP>0), mapping=aes(long,lat,group=group,fill=gCperm2)) +
  geom_path(regions.df3, mapping=aes(long,lat,group=group), col='gray100')+
  scale_fill_gradientn(colors = cols, name = expression(POC~exported~gC~m^-2~y^-1~" "),
                       #breaks=c(100, 200, 300, 400, 500, 600))+
                       breaks=c(100, 200, 300, 400, 500), 
                       labels=c('100', '200', '300', '400', '>500'))+
  scale_x_longitude(breaks = seq(-180,180,45))+
  scale_y_latitude(breaks = seq(-80,80,40))+
  xlab(expression(paste('longitude (',~degree,')',sep='')))+ylab(expression(paste('latitude (',~degree,')',sep='')))+
  #ylim(-90,88)+
  geom_sf(data=World.land, col='black', size=0.1,fill="grey100")+
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major = element_blank(),
        plot.background = element_blank(), legend.position = c(0.5,1.1),legend.direction = "horizontal", plot.margin = unit(c(2,.1,.1,.1), "cm"), 
        legend.margin=margin(t = 0, unit='cm'), legend.title = element_text(size = 12),legend.text = element_text(size = 12), 
        legend.key.width =  unit(2,'cm'))

map.gC

setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")

tiff('Outputs/map.gC.tif', width = 30, height = 15, units = 'cm', res=300)
map.gC
dev.off()


#carbon sequestration using area of ecoregion

regions.df3=regions.df2

quantile(regions.df3$TgC,na.rm = T,.995)
#one region has very large TgC 12, the next largest is 3, smallest is orders of magnitude less.  

#log10 for plotting
regions.df3$TgC=log10(regions.df3$TgC)
regions.df3= regions.df3[!is.na(regions.df3$TgC),]

map.CseqTg=ggplot() +  
  geom_polygon(data=regions.df3%>%filter(regions.df3$NPP>0)%>%filter(regions.df3$area.corrected.P.rock!=0), 
               mapping=aes(long,lat,group=group,fill=TgC)) +
  geom_path(regions.df3, mapping=aes(long,lat,group=group), col='gray100')+
  scale_fill_gradientn(colors = cols, name = expression(TgC~exported~y^"-1 "), #POC~exported~TgC~y^"-1 "
                       breaks=c(-3, -2, -1, 0, 1,2), labels=c(10^-3, 10^-2, 10^-2, 0.1,round(10^1.0,1),round(10^2.0,1)))+
  scale_x_longitude(breaks = seq(-180,180,45))+
  scale_y_latitude(breaks = seq(-80,80,40))+
  xlab(expression(paste('longitude (',~degree,')',sep='')))+ylab(expression(paste('latitude (',~degree,')',sep='')))+
  #ylim(-90,88)+
  geom_sf(data=World.land, col='black', size=0.1,fill="grey100")+
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major = element_blank(),
        plot.background = element_blank(), legend.position = c(0.5,1.1),legend.direction = "horizontal", plot.margin = unit(c(2,.1,.1,.1), "cm"), 
        legend.margin=margin(t = 0, unit='cm'), legend.title = element_text(size = 12),legend.text = element_text(size = 12), 
        legend.key.width =  unit(2,'cm'))
map.CseqTg

setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")

tiff('Outputs/map.CseqTG.tif', width = 30, height = 15, units = 'cm', res=300)
map.CseqTg
dev.off()



# -------------------------------------------------#
#### write main figures of mapped export ####               
# -------------------------------------------------#

figEcoregions=ggarrange(map.export,map.gC,map.CseqTg, nrow = 3, labels = c('A', 'B', 'C'),hjust=-1, font.label = list(face='plain'))

tiff('Outputs/maps.ecoregion.v2.tif', width = 30, height = 45, units = 'cm', res=300)
figEcoregions
dev.off()


figEcoregionsab=ggarrange(map.export,map.gC, nrow = 2, labels = c('A', 'B'),hjust=-1, font.label = list(face='plain'))

tiff('Outputs/maps.ecoregion.percent export and tgC.tif', width = 30, height = 30, units = 'cm', res=300)
figEcoregionsab
dev.off()


# -------------------------------------------------#
#### groundtruthing from regional studies ####               
# -------------------------------------------------#


eco2$export[eco2$marineEco=='West Greenland Shelf']
eco2$exportcombinedSD[eco2$marineEco=='West Greenland Shelf']/sqrt(eco2$N[eco2$marineEco=='West Greenland Shelf'])

#van dheer mheen WA
eco2$export[eco2$marineEco=='Houtman']
eco2$exportcombinedSD[eco2$marineEco=='Houtman']
25*.71

#broch et al. <5%
eco2$export[eco2$marineEco=='Southern Norway']
eco2$exportcombinedSD[eco2$marineEco=='Southern Norway']

colnames(regions.df)
colnames(eco)

#filbee-dexter et al. 2021 Norway model
eco2$export[eco2$marineEco=='Northern Norway and Finnmark']
#percent export
eco2$exportcombinedSD[eco2$marineEco=='Northern Norway and Finnmark']
eco2$exportcombinedSD[eco2$marineEco=='Northern Norway and Finnmark']/
  sqrt(eco2$N[eco2$marineEco=='Northern Norway and Finnmark'])



# -------------------------------------------------#
#### plot by nation ####               
# -------------------------------------------------#

library(maptools)
data(wrld_simpl)

#name of the countries (to make sure the EEZ names match)
wrld_simpl@data$NAME

# setup base map
#get list of country names for each EEZ name
exportbynation<-read_excel('resultsExport.EEZ.merged.xlsx', sheet=4)
#rename North Korea
exportbynation$Country[exportbynation$eez=='North Korean Exclusive Economic Zone']<-'North Korea'
exportbynation<-exportbynation[,3:4]

#merge country names with the export results and make simple data frame for map
exportbynation=merge(df, exportbynation, by='eez')


countries <- exportbynation$Country
metric=exportbynation$TgC
country_df <- data.frame(countries, metric)
country_df$metric=as.numeric(country_df$metric)
country_df=country_df[!is.na(country_df$countries),]

world_map <- map_data(map = "world")

#make sure the USA names align
usa <- map_data('usa')
uk_sf <- ne_states(country = "united kingdom", returnclass = "sf")

#Rename South Korea
country_df$countries[country_df$countries=='Republic of Korea']<-"South Korea"

#plot TgC by nation
g2=ggplot(country_df) +
  geom_map(aes(map_id = countries, fill = metric), map = world_map) +
  
  labs(fill = expression(Exported~POC~(TgC~y^-1)~" ")) +
  geom_polygon(data=usa, aes(x=long, y=lat, group=group), fill='#66C2A5') + 
  geom_sf(data = uk_sf, aes(), fill='#D53E4F',color = NA)+
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), colour = 'black', fill = NA) +
  expand_limits(x = world_map$long, y = world_map$lat) +
  scale_fill_gradientn(colors=brewer.pal(11, "Spectral")[3:11], 
                       breaks=c(0,1,2,3,4,5), labels=c('0','1', '2', '3', '4','5'))+
  theme_void()+
  theme(legend.position = c(0.5,1.1),legend.direction = "horizontal", plot.margin = unit(c(2,.1,.1,.1), "cm"), 
        legend.margin=margin(t = 0, unit='cm'), legend.title = element_text(size = 12),legend.text = element_text(size = 12), 
        legend.key.width =  unit(2,'cm'))


getwd()
dev.off()
tiff('Outputs/maps.Cexport by nation.tif', width = 30, height = 15, units = 'cm', res=300)
g2
dev.off()

