library(sp)
library(raster)
library(rgdal)
library(dismo)
library(dismotools)
Sys.setenvSys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jdk-10.0.2') #change as appropriate
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
library(mistnet)
library(dplyr)
library(gjam)

#Fine scale MN co-occurrence data

load("StateLinesObj.rda")
MN.extent <- extent(-97.5, -89.5, 43.25, 49.55)
MN.state.lines <- crop(states.lines.crop, MN.extent)

hum.jap.MN.extent <- c(-92.6, -91.2, 43.6, 44)

load("HumJap_gjamOutput_relevebased_17Dec2018_setup.rda")

start.time <- Sys.time()
out.hum.jap <-  gjam(~ mean.temp + temp.day.range + ann.precip + precip.cold.q, xdata, ydata, modelList)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

save(xdata, ydata, elist, rlist, modelList, out.hum.jap, file = "HumJap_gjamOutput_relevebased_17Dec2018_gjamOutput.rda")

#Prediction of all of MN
#bioclim.hj <- brick(bioclim.hj)
#vals.new.data <- data.frame(bioclim.hj@data@values)
vals.new.data <- data.frame(values(bioclim.hj))
r <- raster(bioclim.hj, layer = 1)
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

hj.predictions <- as.vector(test.predict$prPresent[,"Humulusjaponicus"])
hj.rast <- raster(nrow=nrow(bioclim.hj), ncol=ncol(bioclim.hj))
extent(hj.rast) <- extent(bioclim.hj)
values(hj.rast) <- hj.predictions
writeRaster(hj.rast, filename = "Hum_Jap_MN_predBox_RedClim__18Dec2018_relevelbased.grd", overwrite=T)

#Plot predictions
v <- viridis(256)
my.at = seq(0,1,1/256)
levelplot(hj.rast, margin = FALSE, maxpixels = 1e5, col.regions=v, at = my.at) + spplot(states.lines.crop, col.regions = "slateblue1")

