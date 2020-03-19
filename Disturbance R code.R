library(rstudioapi)

# Set Working Directory to the directory of this file
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("C:/RobsR/FLORIDA/1 Niche 2020")

library(raster)
library(rgdal)
library(stringr)
library(rworldmap)
library(gbm)
library(SDMTools)
library(stats)
library(spatialEco)
library(corrplot)
library(stats)
library(spatialEco)

source("brt.functions.R")

utm <- "+proj=utm +zone=17 +units=km +datum=WGS84"
utmMeters <- "+proj=utm +zone=17 +units=m +datum=WGS84"

acerRaster <- raster("ACERpredictions.asc")
proj4string(acerRaster) <- utm
acerRasterMeters <- projectRaster(acerRaster, crs=utmMeters)

library(landscapemetrics)

par(mfrow=c(4, 2))

lim <- 20000
maxc <- 0
for (cutoff in seq(0.35, 0.7, 0.05)) {
  acerRasterMetersInt <- acerRasterMeters > cutoff
  enn <- lsm_p_enn(acerRasterMetersInt)
  enn <- subset(enn, class == 1)
  vals <- enn$value
  vals <- vals[vals<=lim]
  histv <- hist(vals, breaks=seq(0, lim, 1000), plot=FALSE)
  maxc <- max(maxc, histv$counts)
}
for (cutoff in seq(0.35, 0.7, 0.05)) {
  acerRasterMetersInt <- acerRasterMeters > cutoff
  enn <- lsm_p_enn(acerRasterMetersInt)
  enn <- subset(enn, class == 1)
  vals <- enn$value
  vals <- vals[vals<=lim]
  hist(vals, main=toString(cutoff), xlab="Distance (m)", ylim=c(0, maxc), breaks=seq(0, lim, 1000))
}

