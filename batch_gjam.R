setwd("c:/users/thomas/desktop/MITPPC_remote_test/MITPPC_remote_test")

library(sp)
library(raster)
library(rgdal)
library(dismo)
library(dismotools)
Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jdk-10.0.2') #change as appropriate
library(rJava)
library(rgeos)
library(viridis)
library(maptools)
library(ggplot2)
library(gbm)
library(alphahull)
library(rasterVis)
library(lattice)
library(rgbif)
library(letsR)
#library(mistnet)
library(dplyr)
library(gjam)
library(parallel)

load("E:/MITPPC Proj/State_Shapes/StateLinesObj.rda")

load("HumJap_gjamOutput_relevebased_17Dec2018_gjamOutput.rda")

bioclim.hj <- brick("bioclim_hj_std.grd")

# read each piece back in R
list2 <- list()
for(i in 1:81){ # change this 9 depending on your number of pieces
  rx <- brick(paste("SplitRas",i,".tif",sep=""))
  list2[[i]] <- rx
}


for(i in 66:70){
  
  vals.new.data <- data.frame(values(list2[[i]]))
  r <- raster(list2[[i]], layer = 1)
  r.pts <- rasterToPoints(r)
  lon.lat <- r.pts[,1:2]
  new.data <- cbind(lon.lat, vals.new.data)
  colnames(new.data) <- c("lon", "lat", "mean.temp", "temp.day.range", "ann.precip", "precip.cold.q")
  n.xdata <- list(xdata = new.data, nsims = 50)
  
  start.time.pred <- Sys.time()
  test.predict <- gjamPredict(out.hum.jap, newdata = n.xdata)
  end.time.pred <- Sys.time()
  time.taken.pred <- end.time.pred - start.time.pred
  time.taken.pred
  
  # # mosaic them, plot mosaic & save output
  # list2$fun   <- max
  # rast.mosaic <- do.call(mosaic,list2)
  # plot(rast.mosaic,axes=F,legend=F,bty="n",box=FALSE)
  # writeRaster(rast.mosaic,filename=paste("Mosaicked_ras",sep=""),
  #             format="GTiff",datatype="FLT4S",overwrite=TRUE)
  
  hj.predictions <- as.vector(test.predict$prPresent[,"Humulusjaponicus"])
  hj.rast <- raster(nrow=nrow(list2[[i]]), ncol=ncol(list2[[i]]))
  extent(hj.rast) <- extent(list2[[i]])
  values(hj.rast) <- hj.predictions
  fname = paste0("Hum_Jap_MN_predBox_",i)
  writeRaster(hj.rast, filename = fname, overwrite=T)
  
}


list3 <- list()
for(i in 1:81){ # change this 9 depending on your number of pieces
  rx <- raster(paste("Hum_Jap_MN_predBox_",i,".grd",sep=""))
  list3[[i]] <- rx
}
# mosaic them, plot mosaic & save output
list3$fun   <- max
rast.mosaic <- do.call(mosaic,list3)
plot(rast.mosaic,axes=F,legend=F,bty="n",box=FALSE)
writeRaster(rast.mosaic,filename=paste("Mosaicked_ras",sep=""),
            format="GTiff",datatype="FLT4S",overwrite=TRUE)


v <- viridis(256)
my.at = seq(0,1,1/256)
levelplot(rast.mosaic, margin = FALSE, maxpixels = 1e5, col.regions=v, at = my.at) + spplot(states.lines.crop, col.regions = "slateblue1")





















