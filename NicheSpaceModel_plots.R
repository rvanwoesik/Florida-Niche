
# Niche Space Model (Plots) -----------------------------------------

library(rstudioapi)

# Set Working Directory to the directory of this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(gbm)
source("brt.functions.R")

windowsFonts(UNI="Arial Unicode MS")

species<-c("ACER", "MYCE", "OFRA", "SSID")
plot.layout <- c(4,3)
cl<-rainbow(length(species))

vars<-c("rangeSST","distFromCoast","noaa_bathymetry","WaveEnergy","Chl.A","MeanSST","Turbidity")
plotLabels<-c("Range of Daily SST (\u2103)","Distance from Coast (km)", "Bathymetry (m)", "Wave Energy (kJ m\u207b\u00b2)", "Chlorophyll-a (mg m\u207b\u00B3)", "Mean Daily SST (\u2103)", "Turbidity (K\u2084\u2089\u2080)(m\u207b\u00B9)")
limits<-data.frame(rbind(rbind(rep(NaN, length(vars)), rep(NaN, length(vars))), rbind(rep(NaN, length(vars)), rep(NaN, length(vars)))), row.names = c("xmin", "xmax", "ymin", "ymax"))
colnames(limits)<-vars

fpredictors <- list(rep(NA, length(species)))
fresponses <- list(rep(NA, length(species)))
fnames <- list(rep(NA, length(species)))

for (s in c(1:length(species))) {
  speciesName<-species[s]
  load(paste0("./results/", speciesName, ".RData"))
  
  pred.names<-resultModel$gbm.call$predictor.names
  dataframe.name <- resultModel$gbm.call$dataframe
  data <- eval(parse(text = dataframe.name))
  nt <- resultModel$n.trees
  predictors <- list(rep(NA,length(pred.names)))
  responses <- list(rep(NA,length(pred.names)))
  
  for (j in c(1:length(pred.names))) {
    k <- match(resultModel$contributions$var[j],pred.names)
    pred.data <- data[,resultModel$gbm.call$gbm.x[k]]
    response.matrix <- plot.gbm(resultModel, i.var = k, n.trees = nt, return.grid = TRUE)
    predictors[[j]] <- response.matrix[,1]
    if (is.factor(data[,resultModel$gbm.call$gbm.x[k]])) {
      predictors[[j]] <- factor(predictors[[j]],levels = levels(data[,resultModel$gbm.call$gbm.x[k]]))
    }
    responses[[j]] <- response.matrix[,2] - mean(response.matrix[,2])
    l<-match(resultModel$contributions$var[j],vars)
    if (is.nan(limits[1,l])) {
      limits[1,l] = min(predictors[[j]])
      limits[2,l] = max(predictors[[j]])
      limits[3,l] = min(responses[[j]])
      limits[4,l] = max(responses[[j]])
    } else {
      limits[1,l] = min(limits[1,l],min(predictors[[j]]))
      limits[2,l] = max(limits[2,l],max(predictors[[j]]))
      limits[3,l] = min(limits[3,l],min(responses[[j]]))
      limits[4,l] = max(limits[4,l],max(responses[[j]]))
    }
  }
  fpredictors[[s]] <- predictors
  fresponses[[s]] <- responses
  fnames[[s]] <- resultModel$contributions$var
}

jpeg(filename="results/top4PD.jpeg", width = 960*4, height = 960*5,
     units = "px", pointsize = 18*4.2, quality = 99, bg = "white")
par(mfrow = plot.layout)
for (j in c(1:length(vars))) {
  plot(0,0, xlab=plotLabels[j], ylab="fitted function",xlim = limits[1:2,j], ylim = limits[3:4,j],type = "n", family="UNI")
  pcl <- c()
  pleg <- c()
  for (s in c(1:length(species))) {
    k <- match(vars[j], fnames[[s]])
    if (is.finite(k)) {
      pcl <- rbind(pcl, s)
      pleg <- rbind(pleg, species[s])
      lines(fpredictors[[s]][[k]],fresponses[[s]][[k]],col = "black",type = 'l', lty=s, lwd=3)
    }
  }
  if (j == 1) {
    legend("topleft", lty = pcl, legend = pleg, lwd=3, cex = 0.8)
  }
}
dev.off()