######################################################################################
# This script estimates the ventilation rates for seaweed POC 
#working directory
setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")
#dont need this
# written by K Filbee-Dexter and A. Pessarrodona  for the EUROMARINE project 21.3.22
######################################################################################

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

#uses years from Siegel et al. 2021
# #get values  ------------------------------------------------------------
setwd("C:/Users/Karen/Dropbox/Euromarine/Global fate of macroalgae 5th paper/Analysis KFD")

#### #get the sequestration
seq <- read_excel("Outputs/sequestration.horizon.ecoregion.xlsx" )
seq$gC<-as.numeric(seq$gC)
seq$export<-as.numeric(seq$export, na.rm = T)
seq$area<-as.numeric(seq$area, na.rm = T)

seq$Cdeep=

totalseq<-sum(seq$gC, na.rm = T)
totalseq


#summary

seq %>% 
  dplyr::select('marineProv', 'Aerobic.anerobic', 'Refractory_Percent') %>% 
  group_by(Tax_group) %>% summarise(count = n(), Mean_refractory=mean(Refractory_Percent, na.rm = T), 
                                    se=sd(Refractory_Percent, na.rm = T)/sqrt(count))

ggplot()
