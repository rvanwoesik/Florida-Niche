
# Niche Space Model (Plots) -----------------------------------------

library(rstudioapi)

# Set Working Directory to the directory of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(raster)
library(rgdal)
library(stringr)
library(gbm)
library(plotKML)
source("brt.functions.R")

longlat<-"+proj=longlat +ellps=WGS84 +no_defs +datum=WGS84"
utm<-"+proj=utm +zone=17 +units=km +datum=WGS84"

FKNMS<-readOGR("./polygons","CoralReef_1000_Dissolve")

FRRPspecies<-read.csv("trainingData/FRRPspecies.csv", header=FALSE)

for (s in 1:nrow(FRRPspecies)) {
  speciesName<-FRRPspecies[[1]][s]
  validSpecies<-as.numeric(unlist(strsplit(as.character(FRRPspecies[[3]][s]), split="/")))
  FRRP<-read.csv("./trainingData/AllCorals2016B.csv")
  FRRP$Year<-as.numeric(str_sub(FRRP$Date,-4,-1))
  FRRP<-subset(FRRP, FRRP$Year>=2011 & FRRP$Year<=2015)
  presence_FRRP<-subset(FRRP, FRRP$Species %in% validSpecies)
  presence_FRRP<-subset(presence_FRRP, !duplicated(presence_FRRP$Code))
  absence_FRRP<-subset(FRRP, !(FRRP$Code %in% presence_FRRP$Code))
  absence_FRRP<-subset(absence_FRRP, !duplicated(absence_FRRP$Code))
  presencePercentage<-(nrow(presence_FRRP)/(nrow(presence_FRRP)+nrow(absence_FRRP)))*100
  if (presencePercentage>10 || speciesName == "ACER") {
    load(paste0("./results/", speciesName, ".RData"))
    
    predictionRaster<-raster(paste0("results/", speciesName, "predictions.asc"))
    proj4string(predictionRaster) <-utm
    longlatPredictions<- projectRaster(predictionRaster, crs = longlat)
    EX <- extent(-83.25548, -79.8523,24.23897, 27.4085)
    extent(longlatPredictions) <- EX
    longlatPredictions <- setExtent(longlatPredictions, EX, keepres=TRUE)
    longlatPredictions<-mask(longlatPredictions, FKNMS)
    
    # as.factor(longlatPredictions@data@names[1])
    # r,raster_name = 'layer.png', colour = r@data@values, file='ms.kml'
    # longlatPredictions@data@values
    kml(longlatPredictions, raster_name = paste0(speciesName, ".png"), file.name = paste0(speciesName, ".kml"), colour = longlatPredictions@data@values)
  }
}
