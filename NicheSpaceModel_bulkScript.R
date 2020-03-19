
# Niche Space Model (FRRP Looped Run) -----------------------------------------

library(rstudioapi)

# Set Working Directory to the directory of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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

windowsFonts(UNI="Arial Unicode MS")



# Initialization ----------------------------------------------------------

# Print And Debug Options
printInfo<-TRUE
debugInfo<-TRUE

# Threshold
#   Defined as prevalence which we force to 0.5 with site weighting.
threshold<-0.5

# Clear Results Folder
unlink("results", recursive=TRUE)

# Create Results Folder
dir.create("results")

# Establish Projections
eqc<-"+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
longlat<-"+proj=longlat +ellps=WGS84 +no_defs"
utm<-"+proj=utm +zone=17 +units=km +datum=WGS84"
aea<-"+proj=aea +lat_1=24 +lat_2=31.5 +lat_0=24 +lon_0=-84 +x_0=400000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"

# Initial Shapefiles
critical<-spTransform(readOGR("./polygons","Acropora_Critical_Habitat"), utm)
subregions<-readOGR("./polygons","Unified_Florida_Reef_Map_v2.0_Regions")
subregionslonglat<- spTransform(subregions, longlat)

# Extent Mask
emptyRaster<-raster(resolution = c(1, 1), xmn = 281.797, xmx =  608.8, ymn= 2689., ymx = 3025.68, crs = utm)
emptyRaster[]<-1
extentmask<-mask(emptyRaster, critical)

# Coral Habitat Extent
FKNMS<-readOGR("./polygons","CoralReef_1000_Dissolve")

# Predictor Data
#(i) depth
#(ii) mean SST
#(iii) variance of SST
#(iv) min SST
#(v) max SST
#(vi) range SST
#(vii) chlorophyll-a concentrations
#(viii) turbidity 
#(ix) wave energy
#(x) distance from coast
#(xi) noaa bathymetry

#(i) depth
depthmr<-raster("./predictorData/depthm.txt", values=TRUE)
proj4string(depthmr)<-aea
depth<-resample(projectRaster(depthmr, crs=utm), extentmask, method="ngb")

#(ii) mean of daily sea surface temperature
meanDailySST<-resample(projectRaster(raster("./predictorData/meanDailySST2011-2015.tif"), crs=utm), extentmask, method="bilinear")

#(iii) variance of daily sea surface temperature
varDailySST<-resample(projectRaster(raster("./predictorData/varianceDailySST2011-2015.tif"), crs=utm), extentmask, method="bilinear")

#(iv) min of daily sea surface temperature
minDailySST<-resample(projectRaster(raster("./predictorData/minDailySST.tif"), crs=utm), extentmask, method="bilinear")

#(v) max of daily sea surface temperature
maxDailySST<-resample(projectRaster(raster("./predictorData/maxDailySST.tif"), crs=utm), extentmask, method="bilinear")

#(vi) range of daily sea surface temperature
rangeDailySST<-resample(projectRaster(raster("./predictorData/rangeDailySST.tif"), crs=utm), extentmask, method="bilinear")
calcRangeDailySST<-maxDailySST-minDailySST
diffRangeDailySST<-rangeDailySST-calcRangeDailySST

#(vii) chlorophyll-a concentrations
meanChla2011_2015<-raster("./predictorData/meanDailyChla2011-2015_1km.asc")
proj4string(meanChla2011_2015)<-eqc
meanChla<-resample(projectRaster(meanChla2011_2015, crs=utm), extentmask, method="bilinear")

#(viii) turbidity
kd490<-resample(projectRaster(raster("./predictorData/avgDailyKD4902013-2015.tif"), crs=utm), extentmask, method = "bilinear")

#(ix) wave energy
waveEnergy<-resample(projectRaster(raster("./predictorData/AverageDailyWE2011-2015.tif"), crs=utm)/1000, extentmask, method="bilinear")

#(x) distance from coast
FLMap<-spTransform(crop(getMap(resolution="high"), c(-86,-77, 23, 29)), longlat)
r2 <- rasterize(spTransform(FLMap, utm), extentmask, 1)
distFromCoast<-raster::distance(r2)

#(xi) noaa bathymetry
noaa_bathy<-resample(projectRaster(raster("./predictorData/noaa_bathymetry.tif"), crs=utm), extentmask, method="bilinear")

# Full Predictor Combination
allPredictors<-stack(meanDailySST, varDailySST,  minDailySST, maxDailySST, rangeDailySST, meanChla, kd490, waveEnergy, distFromCoast, noaa_bathy)
names(allPredictors)<-c("MeanSST", "VarianceSST","minSST","maxSST","rangeSST","Chl-A", "Turbidity","WaveEnergy", "distFromCoast","noaa_bathymetry")
allPredictorsPlotLabels<-c("Mean Daily SST (\u2103)", "Variance of Daily SST (\u2103)", "min()", "max()", "Range of Daily SST (\u2103)", "Chlorophyll-a (mg m\u207b\u00B3)", "Turbidity (K\u2084\u2089\u2080)(m\u207b\u00B9)", "Wave Energy (kJ m\u207b\u00b2)", "Distance from Coast (km)", "Bathymetry (m)")

# Drop useless environment layers
allPredictorsOriginal<-allPredictors
allPredictors<-dropLayer(allPredictors, "minSST")
allPredictors<-dropLayer(allPredictors, "maxSST")
allPredictors<-dropLayer(allPredictors, "VarianceSST")
k <- match(c("minSST", "maxSST", "VarianceSST"), names(allPredictorsOriginal))
allPredictorsPlotLabels<-allPredictorsPlotLabels[-k]

# Save Correlation
par(mar=c(2,2,2,2))
predictor_correlations=layerStats(allPredictors, 'pearson', na.rm=T)
corr_matrix=predictor_correlations$'pearson correlation coefficient'
jpeg(filename="results/predictorCorrelations.jpeg", width = 960*4, height = 960*4,
     units = "px", pointsize = 18*4, quality = 99, bg = "white")
corrplot(corr_matrix, method="color", addCoef.col = "black",sig.level = 0.05)
dev.off()

# Species Solver Loop -------------------------------------------------------------

# RANDOMIZATION
# This program uses randomization. For repeatability, use set.seed(value) to
# establish a consistent randomizer seed.
set.seed(1642320)

# Import Species List Info
FRRPspecies<-read.csv("trainingData/FRRPspecies.csv", header=FALSE)

# Solver Log
logConnection<-NULL
if (debugInfo) {
  logConnection<-file("results/solverLog.txt", open="wt")
}

# nrow(FRRPspecies)
for (s in 1:nrow(FRRPspecies)) {
  
  # Species List (LOOP Head)
  speciesName<-FRRPspecies[[1]][s]
  speciesFullName<-FRRPspecies[[2]][s]
  validSpecies<-as.numeric(unlist(strsplit(as.character(FRRPspecies[[3]][s]), split="/")))
  
  # Solver Log - Species Start
  if (printInfo) {
    cat(toString(speciesFullName), " - Initialization [", toString(s), "]\n", sep="")
  }
  if (debugInfo) {
    writeLines(c("", paste0(speciesName, " - Initialization")), logConnection)
  }
  
  # Full Data Set For Given Species List (FRRP)
  FRRP<-read.csv("./trainingData/AllCorals2016B.csv")
  FRRP$Year<-as.numeric(str_sub(FRRP$Date,-4,-1))
  FRRP<-subset(FRRP, FRRP$Year>=2011 & FRRP$Year<=2015)
  presence_FRRP<-subset(FRRP, FRRP$Species %in% validSpecies)
  presence_FRRP<-subset(presence_FRRP, !duplicated(presence_FRRP$Code))
  absence_FRRP<-subset(FRRP, !(FRRP$Code %in% presence_FRRP$Code))
  absence_FRRP<-subset(absence_FRRP, !duplicated(absence_FRRP$Code))
  presence_Points<-spTransform(SpatialPoints(cbind(presence_FRRP$Longitude, presence_FRRP$Latitude), CRS(longlat)), utm)
  absence_Points<-spTransform(SpatialPoints(cbind(absence_FRRP$Longitude, absence_FRRP$Latitude), CRS(longlat)), utm)
  full_Points<-data.frame(rbind(
    cbind(rep(1, nrow(presence_Points@coords)), presence_Points@coords[,1], presence_Points@coords[,2]),
    cbind(rep(0, nrow(absence_Points@coords)), absence_Points@coords[,1], absence_Points@coords[,2])
  ))
  names(full_Points)<-c("presence", "longitude", "latitude")
  coordinates(full_Points)<-~longitude+latitude
  proj4string(full_Points)<-utm
  
  # Ensure 10% presence sites (or ACER)
  presencePercentage<-(nrow(presence_FRRP)/(nrow(presence_FRRP)+nrow(absence_FRRP)))*100
  if (presencePercentage<10 && speciesName != "ACER" && speciesName != "APAL") {
    # Solver Log - Less Than 10%
    if (debugInfo) {
      writeLines(paste0(speciesName, " - ABORT: Only ", toString(presencePercentage), "% Presence"), logConnection)
    }
    next
  }
  
  # Solver Log - Presence Percentage
  if (debugInfo) {
    if (presencePercentage>=10) {
      writeLines(paste0(speciesName, " - PASS: ", toString(presencePercentage), "% Presence"), logConnection)
    } else {
      writeLines(paste0(speciesName, " - PASS: Species Exception: ", toString(presencePercentage), "% Presence"), logConnection)
    }
  }
  
  # Extract Predictors With Data Points
  allData<-extract(allPredictors, full_Points)
  allData<-as.data.frame(cbind(presence=full_Points$presence,longitude=full_Points@coords[,1], latitude=full_Points@coords[,2], allData))
  allData<-allData[complete.cases(allData),]
  
  # Partition Training and Test Data (RANDOMIZATION)
  sampleSize<-floor(0.8*nrow(allData))
  indices<-sample(seq_len(nrow(allData)),size=sampleSize)
  trainingData<-allData[indices,]
  testData<-allData[-indices,]
  presenceCount<-sum(trainingData$presence==1)
  
  # Ensure Valid Presence In Training Set
  if (presenceCount <= 0) {
    if (debugInfo) {
      writeLines(paste0(speciesName, " - ABORT: No Training Set Presence"), logConnection)
    }
    next
  }
  
  # Ensure Valid Presence In Test Set
  if (sum(testData$presence==1) <= 0) {
    if (debugInfo) {
      writeLines(paste0(speciesName, " - ABORT: No Test Set Presence"), logConnection)
    }
    next
  }
  
  # Establish Site Weights
  weights<-rep(1, nrow(trainingData))
  presenceCount<-sum(trainingData$presence==1)
  weights[trainingData$presence!=1] = rep(presenceCount/(nrow(trainingData)-presenceCount), nrow(trainingData)-presenceCount)
  
  # Prepare For Predictor Loop
  if (printInfo) {
    cat("Initializing Predictor Loop\n")
  }
  
  # Loop Result And Test List
  lpcnt<-0
  result<-(NaN)
  resultValue<-(-Inf)
  resultModel<-NULL
  testList<-list(c(4:length(allData)))
  testListFinalCheck<-list(FALSE)
  usedList<-list()
  
  while (length(testList) > 0) {
    
    # Predictor List (LOOP Head)
    lpcnt<-lpcnt+1
    predictors<-testList[1]
    predictors<-predictors[[1]]
    finalCheck<-testListFinalCheck[1]
    finalCheck<-finalCheck[[1]]
    testList<-testList[-1]
    testListFinalCheck<-testListFinalCheck[-1]
    
    # Optional Output
    if (printInfo) {
      if (finalCheck) {
        cat("Final Testing: [", toString(predictors), "]\n", sep="")
      } else {
        cat("Testing: [", toString(predictors), "]\n", sep="")
      }
    }
    
    {
      sink("nul")
      
      attempt<-0
      trees<-NULL
      
      while (attempt<3 && is.null(trees)) {
        
        # Model Builder (RANDOMIZATION)
        modelBRT <- gbm.step(
          data = trainingData,
          gbm.x = predictors,
          gbm.y = 1,
          family = "bernoulli",
          n.trees = 30,
          tree.complexity = min(c(length(predictors)-1, 6)),
          learning.rate = 0.0015,
          bag.fraction = 0.8,
          site.weights = weights
        )
        
        if (!is.null(modelBRT)) {
          trees<-modelBRT[["n.trees"]]
        }
        attempt<-attempt+1
      }
      
      if (is.null(trees)) {
        sink()
        
        # Solver Log - Model Failure
        if (debugInfo) {
          writeLines(paste0(speciesName, " - ABORT: Model Development Failure"), logConnection)
        }
        next
      }
      
      # Predictions
      testPredictions<-predict.gbm(modelBRT, testData[,predictors], trees)
      testPredictionsLogical<-testPredictions >= threshold
      
      # Determine Overall Performance (Area Under Curve)
      performance<-gbm.roc.area(testData$presence, testPredictionsLogical)
      
      sink()
    }
    
    # If Performance Improved
    if (performance > resultValue) {
      
      # Optional Output
      if (printInfo) {
        cat("Update Top Performer:")
        cat("[", toString(predictors), "] = ", performance, "\n", sep="")
      }
      
      # Establish New Baseline
      result<-predictors
      resultValue<-performance
      resultModel<-modelBRT
      
      # Determine Predictor Ranks And State
      importanceNames<-as.character(modelBRT$contributions$var)
      importance<-modelBRT$contributions$rel.inf
      indices<-order(importance)
      cnt<-floor(length(indices)/2)
      fCheck<-FALSE
      
    } else if (!finalCheck) {
      
      # Determine Predictor Ranks And State
      importanceNames<-as.character(modelBRT$contributions$var)
      importance<-modelBRT$contributions$rel.inf
      indices<-order(importance)
      cnt<-1
      fCheck<-TRUE
      
    } else {
      
      cnt<-0
    }
    
    # Update usedList
    usedList<-c(usedList, list(predictors))
    
    # Add To Options List
    if (cnt > 0) {
      for (i in 1:cnt) {
        id<-indices[i]
        rem<-match(importanceNames[id], names(allPredictors))+3
        npreds<-predictors[predictors != rem]
        if (length(npreds) >= 4 && !(list(npreds) %in% testList) && !(list(npreds) %in% usedList)) {
          
          # Optional Output
          if (printInfo) {
            if (fCheck) {
              cat("Add Final Check: [", toString(npreds), "]\n", sep="")
            } else {
              cat("Add Check: [", toString(npreds), "]\n", sep="")
            }
          }
          
          testList<-c(testList, list(npreds))
          testListFinalCheck<-c(testListFinalCheck, list(fCheck))
        }
      }
      if (cnt > 1) {
        npreds<-predictors
        for (i in 1:cnt) {
          if (length(npreds)<=4) {
            break
          }
          id<-indices[i]
          rem<-match(importanceNames[id], names(allPredictors))+3
          npreds<-npreds[npreds != rem]
          if (!(list(npreds) %in% testList) && !(list(npreds) %in% usedList)) {
            
            # Optional Output
            if (printInfo) {
              cat("Add Final Check: [", toString(npreds), "]\n", sep="")
            }
            
            testList<-c(testList, list(npreds))
            testListFinalCheck<-c(testListFinalCheck, list(TRUE))
          }
        }
      }
    }
    
  }
  
  if (is.null(resultModel)) {
    # Solver Log - Model Failure
    if (debugInfo) {
      writeLines(paste0(speciesName, " - ABORT: Model Development Failure"), logConnection)
    }
    next;
  }
  
  # Solver Log - Model Completion
  if (debugInfo) {
    writeLines(paste0(speciesName, " - Model Complete (AUC: ", toString(resultValue), ")"), logConnection)
  }
  
  # Final Model And Outputs
  resultNames = names(allData)[result];
  predictions<-predict.gbm(resultModel, testData[,result], resultModel[["n.trees"]])
  predictionsLogical<-predictions >= threshold
  performance<-gbm.roc.area(testData$presence, predictionsLogical)
  cmatrix<-confusion.matrix(testData$presence, predictions, threshold=threshold)
  cmatrixAccuracy<-accuracy(testData$presence, predictions, threshold=threshold)
  
  cat("Result: ", performance, "\n", sep="")
  
  # Save Prediction Raster And Plot
  preds<-as.data.frame(allPredictors)[result-3]
  preds<-preds[complete.cases(preds),]
  gbm.predict.grids(
    resultModel,
    preds,
    want.grids = T,
    sp.name = paste0("results/", speciesName, "predictions"),
    pred.vec = rep(-9999, allPredictors@ncols*allPredictors@nrows),
    filepath = "",
    num.col = allPredictors@ncols,
    num.row = allPredictors@nrows,
    xll = 0,
    yll = 0,
    cell.size = 1,
    no.data = -9999
  )
  predictionRaster<-raster(paste0("results/", speciesName, "predictions.asc"))
  proj4string(predictionRaster) <-utm
  longlatPredictions<- projectRaster(predictionRaster, crs = longlat)
  EX <- extent(-83.25548, -79.8523,24.23897, 27.4085)
  extent(longlatPredictions) <- EX
  longlatPredictions <- setExtent(longlatPredictions, EX, keepres=TRUE)
  longlatPredictions<-mask(longlatPredictions, FKNMS)
  writeRaster(longlatPredictions, paste0("results/", speciesName, "predictions.tif"), overwrite=TRUE)
  jpeg(filename=paste0("results/", speciesName, "predictions.jpeg"), width = 960, height = 960, 
       units = "px", pointsize = 12, quality = 95, bg = "white")
  plot(longlatPredictions, legend.args = list(text='Probability of Occurrence', side=4, line=3))
  plot(FLMap, col="gray", add=T)
  dev.off()
  
  # Save Prediction Contributors Plot
  k <- match(resultNames, names(allPredictors))
  labels <- allPredictorsPlotLabels[k]
  if (length(resultNames) <= 6) {
    jpeg(filename=paste0("results/", speciesName, "-predictorContributions.jpeg"), width = 960*8, height = 960*4, 
         units = "px", pointsize = 16*8, quality = 99, bg = "white")
    par(mar=c(5,5,1,1), lwd = 4)
    gbm.plot(resultModel, smooth = TRUE, smooth.lwd = 4, write.title = F,  rug.tick = 0.04, rug.lwd = 4, 
             plot.layout = c(2,3), lwd = 4, x.label = labels, family="UNI")
    dev.off()
  } else {
    jpeg(filename=paste0("results/", speciesName, "-predictorContributions.jpeg"), width = 960*8, height = 960*6, 
         units = "px", pointsize = 16*8, quality = 99, bg = "white")
    par(mar=c(5,5,1,1), lwd = 4)
    gbm.plot(resultModel, smooth = TRUE, smooth.lwd = 4, write.title = F,  rug.tick = 0.04, rug.lwd = 4, 
             plot.layout = c(3,3), lwd = 4, x.label = labels, family="UNI")
    dev.off()
  }
  
  # Suitable Habitat Data
  totalArea <- function(x) {
    return( length(x) ) 
  }
  area <- function(x, t) {
    return( length(x[x >= t]) ) 
  }
  percentArea <- function(x, t) {
    if (length(x) == 0) return(0)
    return( 100 * (length(x[x >= t]) / length(x) ) ) 
  }
  subregionTotalArea <- zonal.stats(x=subregionslonglat, y=longlatPredictions, stat=function(x) {totalArea(x)}, trace=FALSE, plot=FALSE)
  subregionArea <- zonal.stats(x=subregionslonglat, y=longlatPredictions, stat=function(x) {area(x, threshold)}, trace=FALSE, plot=FALSE)
  subregionPercentArea <- zonal.stats(x=subregionslonglat, y=longlatPredictions, stat=function(x) {percentArea(x, threshold)}, trace=FALSE, plot=FALSE)
  subregionResults<-cbind(as.character(subregionslonglat$Region), subregionTotalArea, subregionArea, subregionPercentArea)
  subregionSums<- c("Total", sum(as.numeric(subregionResults[,2]), na.rm = T), sum(as.numeric(subregionResults[,3]), na.rm = T), 100 * (sum(as.numeric(subregionResults[,3]), na.rm = T) / sum(as.numeric(subregionResults[,2]), na.rm = T)))
  subregionResults<-rbind(subregionResults, subregionSums)
  
  # Save
  save(speciesName, speciesFullName, result, resultNames, resultModel, trainingData, testData, cmatrix, cmatrixAccuracy, subregionResults, file=paste0("results/", speciesName, ".RData"))

}

# Solver Log Closure
if (debugInfo) {
  close(logConnection)
}

# TEMPORARY
FRRPspecies<-read.csv("trainingData/FRRPspecies.csv", header=FALSE)
cols<-names(allPredictors)
speciesData<-setNames(data.frame(matrix(ncol = length(cols), nrow = 0)), cols)
pweight<-list()
for (s in 1:nrow(FRRPspecies)) {
  speciesName<-FRRPspecies[[1]][s]
  speciesFullName<-FRRPspecies[[2]][s]
  validSpecies<-as.numeric(unlist(strsplit(as.character(FRRPspecies[[3]][s]), split="/")))
  FRRP<-read.csv("./trainingData/AllCorals2016B.csv")
  FRRP$Year<-as.numeric(str_sub(FRRP$Date,-4,-1))
  FRRP<-subset(FRRP, FRRP$Year>=2011 & FRRP$Year<=2015)
  presence_FRRP<-subset(FRRP, FRRP$Species %in% validSpecies)
  presence_FRRP<-subset(presence_FRRP, !duplicated(presence_FRRP$Code))
  absence_FRRP<-subset(FRRP, !(FRRP$Code %in% presence_FRRP$Code))
  absence_FRRP<-subset(absence_FRRP, !duplicated(absence_FRRP$Code))
  presencePercentage<-(nrow(presence_FRRP)/(nrow(presence_FRRP)+nrow(absence_FRRP)))*100
  pweight[1:length(cols)] <- NA
  if (presencePercentage>10 || speciesName == "ACER") {
    load(paste0("./results/", speciesName, ".RData"))
    for (m in 1:length(resultModel[["contributions"]][["var"]])) {
      k <- match(resultModel[["contributions"]][["var"]][m], names(allPredictors))
      pweight[k]<-resultModel[["contributions"]][["rel.inf"]][m]
    }
  }
  speciesData <- rbind(speciesData, pweight)
}
names(speciesData)<-cols
rownames(speciesData)<-FRRPspecies$V1
write.csv(subset(speciesData, rowSums(is.na(speciesData)) != length(cols)), file="./results/allPredictorsData.csv")



# Compile Full Information Set CSV ------------------------------------------------------------------

# Region Name List
regions <- as.character(subregionslonglat$Region)

FRRPspecies<-read.csv("trainingData/FRRPspecies.csv", header=FALSE)
cols<-c("TotalPresence", "TotalAbsence", "TrainingPresence", "TrainingAbsence", "TestPresence", "TestAbsence", "ConfusionP0O0", "ConfusionP1O0", "ConfusionP0O1", "ConfusionP1O1", "AUC", "Sensitivity", "Specificity")
cols<-c(cols, paste0("Area of ", as.character(subregionslonglat$Region)), "Total Area")
cols<-c(cols, paste0("Percent Area of ", as.character(subregionslonglat$Region)), "Total Percent Area")
cols<-gsub(" ", "", cols, fixed = TRUE)
speciesData<-setNames(data.frame(matrix(ncol = length(cols), nrow = 0)), cols)
for (s in 1:nrow(FRRPspecies)) {
  speciesName<-FRRPspecies[[1]][s]
  speciesFullName<-FRRPspecies[[2]][s]
  validSpecies<-as.numeric(unlist(strsplit(as.character(FRRPspecies[[3]][s]), split="/")))

  FRRP<-read.csv("./trainingData/AllCorals2016B.csv")
  FRRP$Year<-as.numeric(str_sub(FRRP$Date,-4,-1))
  FRRP<-subset(FRRP, FRRP$Year>=2011 & FRRP$Year<=2015)
  presence_FRRP<-subset(FRRP, FRRP$Species %in% validSpecies)
  presence_FRRP<-subset(presence_FRRP, !duplicated(presence_FRRP$Code))
  absence_FRRP<-subset(FRRP, !(FRRP$Code %in% presence_FRRP$Code))
  absence_FRRP<-subset(absence_FRRP, !duplicated(absence_FRRP$Code))
  
  trainP <- NA
  trainA <- NA
  testP <- NA
  testA <- NA
  cmatrix <- matrix(data=NA,nrow=2,ncol=2)
  cmatrixAccuracy <- matrix(data=NA,nrow=1,ncol=7)
  areas <- matrix(data=NA,nrow=1,ncol=13)
  areasP <- matrix(data=NA,nrow=1,ncol=13)
  
  presencePercentage<-(nrow(presence_FRRP)/(nrow(presence_FRRP)+nrow(absence_FRRP)))*100
  if (presencePercentage>10 || speciesName == "ACER") {
    load(paste0("./results/", speciesName, ".RData"))
    trainP <- nrow(subset(trainingData, presence == 1))
    trainA <- nrow(subset(trainingData, presence == 0))
    testP <- nrow(subset(testData, presence == 1))
    testA <- nrow(subset(testData, presence == 0))
    areas<-as.numeric(subregionResults[,3])
    areasP<-as.numeric(subregionResults[,4])
  }
  iframe <- data.frame(TotalPresence=nrow(presence_FRRP), TotalAbsence=nrow(absence_FRRP), TrainingPresence=trainP, TrainingAbsence=trainA, TestPresence=testP, TestAbsence=testA, ConfusionP0O0=cmatrix[1,1], ConfusionP1O0=cmatrix[2,1], ConfusionP0O1=cmatrix[1,2], ConfusionP1O1=cmatrix[2,2], AUC=cmatrixAccuracy[2], Sensitivity=cmatrixAccuracy[4], Specificity=cmatrixAccuracy[5])
  for (ss in 1:length(regions)) {
    iframe[paste0("Areaof", as.character(ss))] <- areas[ss]
  }
  iframe["TotalArea"] <- areas[13]
  for (ss in 1:length(regions)) {
    iframe[paste0("PercentAreaof", as.character(ss))] <- areasP[ss]
  }
  iframe["TotalPercentArea"] <- areasP[13]
  speciesData<-rbind(speciesData, setNames(iframe, cols))
}
rownames(speciesData)<-FRRPspecies$V1
write.csv(subset(speciesData, !is.na(speciesData$TrainingAbsence)), file="./results/allData.csv")
