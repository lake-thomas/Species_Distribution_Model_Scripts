#Library List
library(sp)
library(raster)
library(rgdal)
library(dismo)
#library(dismotools)
#dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_161.jdk/Contents/Home/jre/lib/server/libjvm.dylib') #on office computer
#dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_141.jdk/Contents/Home/jre/lib/server/libjvm.dylib') #on home laptop
#Sys.setenv(JAVA_HOME='/Library/Java/JavaVirtualMachines/jdk-11.0.1.jdk')
#library(rJava)
library(rgeos)
library(viridis)
library(maptools)
library(maps)
library(mapdata)
library(httr)
library(ggplot2)
library(prism)
library(gbm)
library(vegan)
library(alphahull)
library(ggplot2)
library(rmaxent)
library(rasterVis)
library(ROCR)
library(biomod2)
library(vcd)
library(ztable)
library(ggpubr) #masks extract and rotate from 'raster'
#library(multipanelfigure)
library(grid)
library(gridExtra)
library(rasterVis)
library(lattice)
library(devtools)
#library(dismotools)
data(wrld_simpl)
data("stateMapEnv")

#Background files for processing and for generating maps

#MN state polygon for mask clipping for pictures
states <- readOGR("D:/MITPPC Proj/State_Shapes/cb_2016_us_state_5m.shp")
MN.state.poly.nad83 <- states[13,]
pred.proj <- '+proj=longlat +datum=WGS84'
MN.state.poly <- spTransform(MN.state.poly.nad83, CRS(pred.proj))

#State lines object for ploting state lines
load("D:/MITPPC Proj/State_Shapes/StateLinesObj.rda")

#Parameters for levelplots

my.at <- seq(0,1,1/256)
my.at.diff <- seq(-1,1,1/128)

p.strip <- list(cex=0.5, lines = 1, col = "gray6")

v <- viridis(256)
m <- magma(256)
pl <- plasma(256)

v.2 = v[129:256]
mw.at = seq(0.5,1,1/128) 

#Extent object for scale of plotting
MN.extent <- extent(-97.5, -89.5, 43.25, 49.55)
MN.state.lines <- crop(states.lines.crop, MN.extent)

#North American Extent
NA.coord.extent <- extent(-130, -50, 8.8, 60)


climatestack <- brick("D:/MITPPC Proj/Climate/NA_wc2_BioClim19.grd")
names(climatestack) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

#Identify each raster layer as an object so that they can be used and manipulated individually

mean.temp.bio1 <- raster(climatestack, layer = 1) #(Celcius)
mean.temp.day.range.bio2 <- raster(climatestack, layer = 2) #Mean diurnal temperature range (mean(period max-min)) (Celcius)
isotherm.bio3 <- raster(climatestack, layer = 3) #Isothermality (Bio2/Bio7)
temp.seas.bio4 <- raster(climatestack, layer = 4) #Temperature seasonality (C of V)
max.temp.warm.w.bio5 <- raster(climatestack, layer = 5) #Maximum temperature of the warmest week (Celcius)
min.temp.cold.w.bio6 <- raster(climatestack, layer = 6) #Minimum temperature of the coldest week (Celcius)
ann.temp.range.bio7 <- raster(climatestack, layer = 7) #Annual temperature range (Bio5-Bio6) (Celcius)
mean.temp.wet.q.bio8 <-raster(climatestack, layer = 8) #Mean temperature of wettest quarter (Celcius)
mean.temp.dry.q.bio9 <- raster(climatestack, layer = 9) #Mean temperature of driest quarter (Celcius)
mean.temp.warm.q.bio10 <- raster(climatestack, layer = 10) #Mean temperature of warmest quarter (Celcius)
mean.temp.cold.q.bio11 <- raster(climatestack, layer = 11) #Mean temperature of the coldest quarter (Celcius)
ann.precip.bio12 <- raster(climatestack, layer = 12) #Annual precipitation (mm)
precip.wet.w.bio13 <- raster(climatestack, layer = 13) #Precipitation of the wettest week (mm)
precip.dry.w.bio14 <-raster(climatestack, layer = 14) #Preciptitation of the driest week (mm)
precip.seas.bio15 <- raster(climatestack, layer = 15) #Preciptitation seasonality
precip.wet.q.bio16 <- raster(climatestack, layer = 16) #Precipitation of the wettest quarter (mm)
precip.dry.q.bio17 <- raster(climatestack, layer = 17) #Preciptitation of the driest quarter (mm)
precip.warm.q.bio18 <- raster(climatestack, layer = 18) #Precipitation of the warmest quarter (mm)
precip.cold.q.bio19 <- raster(climatestack, layer = 19) #Precipitation of the coldest quarter (mm)

#Restack only layers of interest in a new object (Here we exclude the PC variables)

bioclim.only <- stack(mean.temp.bio1, mean.temp.day.range.bio2, isotherm.bio3, temp.seas.bio4, max.temp.warm.w.bio5, min.temp.cold.w.bio6, ann.temp.range.bio7, mean.temp.wet.q.bio8, mean.temp.dry.q.bio9, mean.temp.warm.q.bio10, mean.temp.cold.q.bio11, ann.precip.bio12, precip.wet.w.bio13, precip.dry.w.bio14,precip.seas.bio15, precip.wet.q.bio16, precip.dry.q.bio17, precip.warm.q.bio18, precip.cold.q.bio19)




for(j in 1:3){

setwd("F:/ENMEval Models")

base <- "F:/ENMEval Models/"
folders <- list.files()
fp <- paste(base, folders[j], sep="")
rda_files <- list.files(fp, full.names = TRUE, pattern ="*.rda")

for(i in 1:length(rda_files)){
  
  rda_file <- load(rda_files[i])
  prediction <- rmaxent::project(best.model, bioclim.nb.3)
  
  plot1 <- levelplot(prediction[[3]], margin = FALSE, maxpixels = 1e6, col.regions = v, at = my.at, xlab = list("Longitude", cex = 0.75), ylab = list("Latitude", cex = 0.75), scales = list(cex = 0.5, cex.axis = 0.5, cex.lab = 0.5)) + spplot(states.lines.crop, col.regions = "slateblue1")
  
  name <- paste(substr(rda_files[i], 1, nchar(rda_files[i])-4), ".pdf", sep="")
  
  {cairo_pdf(file = name, width = 8, height = 7, family = "Arial")
    print(plot1)
  }
  dev.off()
  
  
  rm(prediction, best.model, plot1)
  gc()
  
}

}



































