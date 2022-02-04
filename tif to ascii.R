#Libraries
library(gjam)
library(sp)
library(raster)
library(rgdal)
library(dismo)
library(rgeos)
library(viridis)
library(maptools)
library(rasterVis)


#################
###Bittercress###
#################

#Bittercress file path
carimp.jsdm.file.path <- "c:/users/thomas/desktop"

#Bittercress Climate Data for prediction
new.data <- read.csv(paste(carimp.jsdm.file.path, "CarImp_MN_ClimData.csv", sep = "/"))
n.xdata <- list(xdata = new.data, nsims = 50)

#Spp Box JSDM
#load(paste(carimp.jsdm.file.path, "CarImp_gjamOutput_SppBox_19Feb2019.rda", sep = "/"))
#test.predict <- gjamPredict(out.car.imp, newdata = n.xdata)
#save(test.predict, file = paste(carimp.jsdm.file.path, "CarImp_SppBox_gjamPredictObj.rda", sep = "/"))

#Part State JSDM
load(paste(carimp.jsdm.file.path, "CarImp_gjamOutput_PartState_19Feb2019.rda", sep = "/"))
test.predict <- gjamPredict(out.car.imp, newdata = n.xdata)
save(test.predict, file = paste(carimp.jsdm.file.path, "CarImp_PartState_gjamPredictObj.rda", sep = "/"))

#Full State JSDM
load(paste(carimp.jsdm.file.path, "CarImp_gjamOutput_FullState_19Feb2019.rda", sep = "/"))
test.predict <- gjamPredict(out.car.imp, newdata = n.xdata)
save(test.predict, file = paste(carimp.jsdm.file.path, "CarImp_FullState_gjamPredictObj.rda", sep = "/"))


#################
###Bittersweet###
#################

#Bittersweet file path
celorb.jsdm.file.path <- "c:/users/thomas/desktop"

#Bittersweet Climate Data for prediction
new.data <- read.csv(paste(celorb.jsdm.file.path, "CelOrb_MN_ClimData.csv", sep = "/"))
n.xdata <- list(xdata = new.data, nsims = 50)

#Spp Box JSDM
load(paste(celorb.jsdm.file.path, "CelOrb_gjamOutput_SppBox_19Feb2019.rda", sep = "/"))
test.predict <- gjamPredict(out.cel.orb, newdata = n.xdata)
save(test.predict, file = paste(celorb.jsdm.file.path, "CelOrb_SppBox_gjamPredictObj.rda", sep = "/"))

#Part State JSDM
load(paste(celorb.jsdm.file.path, "CelOrb_gjamOutput_PartState_19Feb2019.rda", sep = "/"))
test.predict <- gjamPredict(out.cel.orb, newdata = n.xdata)
save(test.predict, file = paste(celorb.jsdm.file.path, "CelOrb_PartState_gjamPredictObj.rda", sep = "/"))

#Full State JSDM
load(paste(celorb.jsdm.file.path, "CelOrb_gjamOutput_FullState_19Feb2019.rda", sep = "/"))
test.predict <- gjamPredict(out.cel.orb, newdata = n.xdata)
save(test.predict, file = paste(celorb.jsdm.file.path, "CelOrb_FullState_gjamPredictObj.rda", sep = "/"))


##########
###Hops###
##########

#Hops File Path
humjap.jsdm.file.path <- "c:/users/thomas/desktop"

#Hops Climate Data fro prediction
new.data <- read.csv(paste(humjap.jsdm.file.path, "HumJap_MN_ClimData.csv", sep = "/"))
n.xdata <- list(xdata = new.data, nsims = 50)

#Spp Box JSDM
load(paste(humjap.jsdm.file.path, "HumJap_gjamOutput_SppBox_19Feb2019.rda", sep = "/"))
test.predict <- gjamPredict(out.hum.jap, newdata = n.xdata)
save(test.predict, file = paste(humjap.jsdm.file.path, "HumJap_SppBox_gjamPredictObj.rda", sep = "/"))

#Part State JSDM
load(paste(humjap.jsdm.file.path, "HumJap_gjamOutput_PartState_19Feb2019.rda", sep = "/"))
test.predict <- gjamPredict(out.hum.jap, newdata = n.xdata)
save(test.predict, file = paste(humjap.jsdm.file.path, "HumJap_PartState_gjamPredictObj.rda", sep = "/"))

#Full State JSDM
load(paste(humjap.jsdm.file.path, "HumJap_gjamOutput_FullState_19Feb2019.rda", sep = "/"))
test.predict <- gjamPredict(out.hum.jap, newdata = n.xdata)
save(test.predict, file = paste(humjap.jsdm.file.path, "HumJap_FullState_gjamPredictObj.rda", sep = "/"))




library(sp)
library(raster)
library(rgdal)
library(dismo)
library(dismotools)
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
library(rmaxent)
library(rasterVis)
library(ROCR)
library(biomod2)
library(vcd)
library(ztable)
library(ggpubr)
library(grid)
library(lattice)
library(devtools)
#data(wrld_simpl)
#data("stateMapEnv")
library(caret)
Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jdk-10.0.2')
library(rJava)
#library(ENMeval)
library(doParallel)
library(parallel)
library(enmSdm)
library(maxnet)
library(gridExtra)

source('c:/users/thomas/desktop/enmeval_functions.R')

load("D:/MITPPC Proj/State_Shapes/StateLinesObj.rda")



#Converting all .rda files to .asc and .tif files




rda_files <- list.files("F:/ENMEval Models/bittercress", pattern = "*50.tif", full.names = TRUE)

for(i in 1:length(rda_files)){
  
  r <- raster(rda_files[i])
  
  ra <- aggregate(r, fact=2)
  
  filename <- paste(tools::file_path_sans_ext(rda_files[i]), ".asc", sep="")
  
  writeRaster(raster_obj, filename = filename, format="ascii", datatype="INT4S")
  
  print(sprintf("Writing Filename %s", filename))
  
  
}


























































