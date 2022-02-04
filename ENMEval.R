##########################
# Master ENMEvaluate SDM #
##########################


#
#Thomas Lake
#March 29 2019
#

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


###Load All Climate Data

###########
##Bioclim##
###########

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

#Future Climates

###########
###CCSM4###
###########

#read .grd from file
climate.stack.ccsm4_50 <- brick("D:/MITPPC Proj/Future Climates/North America/NA_cc45bi50_BioClim19.grd")
names(climate.stack.ccsm4_50) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")

climate.stack.ccsm4_70 <- brick("D:/MITPPC Proj/Future Climates/North America/NA_cc45bi70_BioClim19.grd")
names(climate.stack.ccsm4_70) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")

##########
###GFDL###
##########

#read .grd from file
climate.stack.gfdl_50 <- brick("D:/MITPPC Proj/Future Climates/North America/NA_gd45bi50_BioClim19.grd")
names(climate.stack.gfdl_50) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")

climate.stack.gfdl_70 <- brick("D:/MITPPC Proj/Future Climates/North America/NA_gd45bi70_BioClim19.grd")
names(climate.stack.gfdl_70) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")

##########
###IPSL###
##########

#read .grd from file
climate.stack.ipsl_50 <- brick("D:/MITPPC Proj/Future Climates/North America/NA_ip45bi50_BioClim19.grd")
names(climate.stack.ipsl_50) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")

climate.stack.ipsl_70 <- brick("D:/MITPPC Proj/Future Climates/North America/NA_ip45bi70_BioClim19.grd")
names(climate.stack.ipsl_70) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")

###############
###MIROC-ESM###
###############

#read .grd from file
climate.stack.miroc_50 <- brick("D:/MITPPC Proj/Future Climates/North America/NA_mr45bi50_BioClim19.grd")
names(climate.stack.miroc_50) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")

climate.stack.miroc_70 <- brick("D:/MITPPC Proj/Future Climates/North America/NA_mr45bi70_BioClim19.grd")
names(climate.stack.miroc_70) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")


###############
###MRI-CGCM3###
###############

#read .grd from file
climate.stack.mri_50 <- brick("D:/MITPPC Proj/Future Climates/North America/NA_mg45bi50_BioClim19.grd")
names(climate.stack.mri_50) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")

climate.stack.mri_70 <- brick("D:/MITPPC Proj/Future Climates/North America/NA_mg45bi70_BioClim19.grd")
names(climate.stack.mri_70) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")



#################################
###Climate Datasets for Species##
#################################

#Common Teasel
bioclim.ct.3 <- stack(temp.seas.bio4, mean.temp.warm.q.bio10, mean.temp.cold.q.bio11, precip.seas.bio15, precip.wet.q.bio16)

#Black Swallowwort
bioclim.sw.3 <- stack(temp.seas.bio4, mean.temp.warm.q.bio10, ann.precip.bio12, precip.seas.bio15, precip.wet.q.bio16)

#Brown Knapweed
bioclim.bk.3 <- stack(mean.temp.bio1, temp.seas.bio4, mean.temp.warm.q.bio10, ann.precip.bio12, precip.seas.bio15, precip.wet.q.bio16)

#Dalmatian Toadflax
bioclim.tf.3 <- stack(temp.seas.bio4, mean.temp.warm.q.bio10, mean.temp.cold.q.bio11, precip.seas.bio15)

#wild parsnip
bioclim.parsnip.3 <- stack(mean.temp.warm.q.bio10, mean.temp.cold.q.bio11, ann.precip.bio12, precip.seas.bio15, precip.wet.q.bio16)

#Common Tansy
bioclim.tansy.3 <- stack(mean.temp.warm.q.bio10, mean.temp.cold.q.bio11, precip.seas.bio15, precip.wet.q.bio16)

#Japanese Hops
bioclim.jh.3 <- stack(temp.seas.bio4, mean.temp.warm.q.bio10, precip.seas.bio15, precip.wet.q.bio16)

#Narrowleaf bittercress
bioclim.nb.3 <- stack(temp.seas.bio4, mean.temp.warm.q.bio10, precip.seas.bio15, precip.wet.q.bio16)

#Oriental Bittersweet
bioclim.ob.3 <- stack(mean.temp.warm.q.bio10, mean.temp.cold.q.bio11, precip.seas.bio15, precip.wet.q.bio16)



env.predictors <- c(bioclim.ct.3, bioclim.sw.3, bioclim.bk.3, bioclim.tf.3, bioclim.parsnip.3, bioclim.tansy.3, bioclim.jh.3, bioclim.nb.3, bioclim.ob.3)


#Load Species Data
#load("E:/MITPPC Proj/State_Shapes/StateLinesObj.rda") #occupies a lot of ram

#read species data
species.data.path <- "D:/MITPPC Proj/Biomod2 Modelling/occurrence data/Final Species Coords/"
ct.file.path <- paste(species.data.path, "common_teasel_coords.csv", sep = "/") #common teasel
sw.file.path <- paste(species.data.path, "swallowwort_dataset_2.csv", sep = "/") #black swallowwort
bk.file.path <- paste(species.data.path, "knapweed_occurrence_data_2018.csv", sep = "/") #brown knapweed
tf.file.path <- paste(species.data.path, "dalmatian_toadflax_coords_counties.csv", sep = "/") #dalmatian toadflax
wp.file.path <- paste(species.data.path, "wild_parsnip_1km.csv", sep = "/")
ctansy.file.path <- paste(species.data.path, "common_tansy_1km.csv", sep = "/")
jh.file.path <- paste(species.data.path, "humulus_japonicus_clean_22June2018.csv", sep = "/")
nb.file.path <- paste(species.data.path, "narrowleaf_bittercress_dataset_28June2018.csv", sep = "/")
ob.file.path <- paste(species.data.path, "oriental_bittersweet_dataset_16May2018.csv", sep = "/")

#model batch parameters

species.names <- c("teasel", "swallowwort", "knapweed", "toadflax", "parsnip", "tansy", "hops", "bittercress", "bittersweet")
gdk.values <- c(0, 1, 3)
downsample.values <- c(0.01, 0.1)
downsample.units <- c("1km", "10km")

#Absolute file paths for species GDK files
teasel_bias_list <- c('', 'c:/users/thomas/desktop/GDK_Master/teasel_gdk_1.tif', 'c:/users/thomas/desktop/GDK_Master/teasel_gdk_3.tif')
swallowwort_bias_list <- c('', 'c:/users/thomas/desktop/GDK_Master/swallowwort_gdk_1.tif', 'c:/users/thomas/desktop/GDK_Master/swallowwort_gdk_3.tif')
knapweed_bias_list <- c('', 'c:/users/thomas/desktop/GDK_Master/knapweed_gdk_1.tif', 'c:/users/thomas/desktop/GDK_Master/knapweed_gdk_3.tif')
toadflax_bias_list <- c('', 'c:/users/thomas/desktop/GDK_Master/toadflax_gdk_1.tif', 'c:/users/thomas/desktop/GDK_Master/toadflax_gdk_3.tif')
parsnip_bias_list <- c('', 'c:/users/thomas/desktop/GDK_Master/parsnip_gdk_1.tif', 'c:/users/thomas/desktop/GDK_Master/parsnip_gdk_3.tif')
tansy_bias_list <- c('', 'c:/users/thomas/desktop/GDK_Master/tansy_gdk_1.tif', 'c:/users/thomas/desktop/GDK_Master/tansy_gdk_3.tif')
hops_bias_list <- c('', 'c:/users/thomas/desktop/GDK_Master/japanese_hops_15_gdk.tif', 'c:/users/thomas/desktop/GDK_Master/japanese_hops_25_gdk.tif')
bittercress_bias_list <- c('', 'c:/users/thomas/desktop/GDK_Master/narrowleaf_bittercress_15_gdk.tif', 'c:/users/thomas/desktop/GDK_Master/narrowleaf_bittercress_25_gdk.tif')
bittersweet_bias_list <- c('', 'c:/users/thomas/desktop/GDK_Master/oriental_bittersweet_15_gdk.tif', 'c:/users/thomas/desktop/GDK_Master/oriental_bittersweet_25_gdk.tif')

bias.list <- c(teasel_bias_list, swallowwort_bias_list, knapweed_bias_list, toadflax_bias_list, parsnip_bias_list, tansy_bias_list, hops_bias_list, bittercress_bias_list, bittersweet_bias_list)

#Absolute file paths for species presence data
species.file.paths <- c(ct.file.path, sw.file.path, bk.file.path, tf.file.path, wp.file.path, ctansy.file.path, jh.file.path, nb.file.path, ob.file.path)





###################
## MASTER RUN #####
###################

#
master_results_dataframe <- data.frame()


for(i in 1:length(species.file.paths)){ #models for all species

  cat(sprintf("------- Starting Species %s --------\n", species.names[i]))

  
  for(j in 1:length(downsample.values)){
  
    cat(sprintf("downsample value %s \n",downsample.values[[j]]))
    
  
    for(k in 1:length(gdk.values)){
    
      cat(sprintf("gdk value %s \n",gdk.values[[k]]))
      
      bias_list <- get(paste(species.names[i], "_bias_list", sep="")) #get the species_bias_list object 
    
      d <- prepareData(species.file.paths[i], downsample_value = downsample.values[j], gdk_value = gdk.values[k], bias_file = bias_list[k])
    
      o <- runENMModel(d[[1]], env.predictors[[i]], d[[2]], species_name = species.names[i], gdk_value = gdk.values[k], downsample_value = downsample.units[j])
      
      master_results_dataframe <- rbind( o, master_results_dataframe)
      
      print(master_results_dataframe)
    
    }
    
  }
  
}

write.csv(master_results_dataframe, file = "F:/master_species_results_dataframe.csv")










#Testing getting species blocks and checkerboards points
#Question: which data is in which blocK?
#Question: do the blocks appear differnet from different background sampling or downsampling values?

#swallowwort

bias_list <- get(paste(species.names[2], "_bias_list", sep="")) #get the species_bias_list object 

#prepare species occurrence and background data
d <- prepareData(species.file.paths[2], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

blocks <- get.block(d[[1]], d[[2]])

occ_blocks <- cbind(d[[1]], blocks$occ.grp)

block1 <- occ_blocks[which(occ_blocks[,3]==1),]

block2 <- occ_blocks[which(occ_blocks[,3]==2),]

block3 <- occ_blocks[which(occ_blocks[,3]==3),]

block4 <- occ_blocks[which(occ_blocks[,3]==4),]

{cairo_pdf(file = "c:/users/thomas/desktop/Black Swallowwort Block Validation.pdf", width = 8, height = 7)
  
  plot(mean.temp.bio1, col="gray90", legend=FALSE, main="Black Swallowwort Block Validation", xlab = "Longitude", ylab = "Latitude")

  points(block4, col="#440154FF", pch=20, cex=0.8)
  points(block3, col="#31688EFF", pch=20, cex = 0.8)
  points(block2, col="#35B779FF", pch=20, cex = 0.8)
  points(block1, col="#FDC726", pch=20, cex = 0.8)
  abline(v=-72)
  abline(h=42.3)
}
dev.off()


#checkerboard

check <- get.checkerboard1(d[[1]], mean.temp.bio1, d[[2]], aggregation.factor = 500)

chk_blocks <- cbind(d[[1]], check$occ.grp)

block1 <- chk_blocks[which(chk_blocks[,3]==1),]

block2 <- chk_blocks[which(chk_blocks[,3]==2),]

block3 <- chk_blocks[which(chk_blocks[,3]==3),]

block4 <- chk_blocks[which(chk_blocks[,3]==4),]

{cairo_pdf(file = "c:/users/thomas/desktop/Black Swallowwort Ckeckerboard Validation.pdf", width = 8, height = 7)
  
  plot(mean.temp.bio1, col="gray90", legend=FALSE, main="Black Swallowwort Checkerboard Validation", xlab = "Longitude", ylab = "Latitude")
  
  points(block2, col="#35B779FF", pch=20, cex = 0.8)
  points(block1, col="#440154FF", pch=20, cex=0.8)
  
  abline(v=-97)
  abline(v=-93)
  abline(v=-89)
  abline(v=-85)
  abline(v=-81)
  abline(v=-77)
  abline(v=-73)
  abline(v=-69)
  
  
  abline(h=35)
  abline(h=39)
  abline(h=43)
  abline(h=47)

}
dev.off()


#teasel


bias_list <- get(paste(species.names[1], "_bias_list", sep="")) #get the species_bias_list object 

#prepare species occurrence and background data
d <- prepareData(species.file.paths[1], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

blocks <- get.block(d[[1]], d[[2]])

occ_blocks <- cbind(d[[1]], blocks$occ.grp)

block1 <- occ_blocks[which(occ_blocks[,3]==1),]

block2 <- occ_blocks[which(occ_blocks[,3]==2),]

block3 <- occ_blocks[which(occ_blocks[,3]==3),]

block4 <- occ_blocks[which(occ_blocks[,3]==4),]

#{cairo_pdf(file = "c:/users/thomas/desktop/Common Teasel Block Validation.pdf", width = 8, height = 7)
  
plot(mean.temp.bio1, col="gray90", legend=FALSE, main="Common Teasel Block Validation", xlab = "Longitude", ylab = "Latitude")
  
points(block4, col="#440154FF", pch=20, cex=0.8)
points(block3, col="#31688EFF", pch=20, cex = 0.8)
points(block2, col="#35B779FF", pch=20, cex = 0.8)
points(block1, col="#FDC726", pch=20, cex = 0.8)
  
#}
#dev.off()


#checkerboard

check <- get.checkerboard1(d[[1]], mean.temp.bio1, d[[2]], aggregation.factor = 500)

chk_blocks <- cbind(d[[1]], check$occ.grp)

block1 <- chk_blocks[which(chk_blocks[,3]==1),]

block2 <- chk_blocks[which(chk_blocks[,3]==2),]

block3 <- chk_blocks[which(chk_blocks[,3]==3),]

block4 <- chk_blocks[which(chk_blocks[,3]==4),]

{cairo_pdf(file = "c:/users/thomas/desktop/Common Teasel Ckeckerboard Validation.pdf", width = 8, height = 7)
  
  plot(mean.temp.bio1, col="gray90", legend=FALSE, main="Black Swallowwort Checkerboard Validation", xlab = "Longitude", ylab = "Latitude")
  
  points(block2, col="#35B779FF", pch=20, cex = 0.8)
  points(block1, col="#440154FF", pch=20, cex=0.8)
  
  abline(v=-97)
  abline(v=-93)
  abline(v=-89)
  abline(v=-85)
  abline(v=-81)
  abline(v=-77)
  abline(v=-73)
  abline(v=-69)
  
  
  abline(h=35)
  abline(h=39)
  abline(h=43)
  abline(h=47)
  
}
dev.off()




#knapweed

bias_list <- get(paste(species.names[3], "_bias_list", sep="")) #get the species_bias_list object 

#prepare species occurrence and background data
d <- prepareData(species.file.paths[3], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

blocks <- get.block(d[[1]], d[[2]])

occ_blocks <- cbind(d[[1]], blocks$occ.grp)

block1 <- occ_blocks[which(occ_blocks[,3]==1),]

block2 <- occ_blocks[which(occ_blocks[,3]==2),]

block3 <- occ_blocks[which(occ_blocks[,3]==3),]

block4 <- occ_blocks[which(occ_blocks[,3]==4),]

plot(mean.temp.bio1, col="gray", legend=FALSE, main="knapweed Block")

points(block4, col="grey", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))

#checkerboard

check <- get.checkerboard2(d[[1]], mean.temp.bio1, d[[2]], 5)

chk_blocks <- cbind(d[[1]], check$occ.grp)

block1 <- chk_blocks[which(chk_blocks[,3]==1),]

block2 <- chk_blocks[which(chk_blocks[,3]==2),]

block3 <- chk_blocks[which(chk_blocks[,3]==3),]

block4 <- chk_blocks[which(chk_blocks[,3]==4),]

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Knapweed Checkerboard")

points(block4, col="grey", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))






#toadflax

bias_list <- get(paste(species.names[4], "_bias_list", sep="")) #get the species_bias_list object 

#prepare species occurrence and background data
d <- prepareData(species.file.paths[4], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

blocks <- get.block(d[[1]], d[[2]])

occ_blocks <- cbind(d[[1]], blocks$occ.grp)

block1 <- occ_blocks[which(occ_blocks[,3]==1),]

block2 <- occ_blocks[which(occ_blocks[,3]==2),]

block3 <- occ_blocks[which(occ_blocks[,3]==3),]

block4 <- occ_blocks[which(occ_blocks[,3]==4),]

plot(mean.temp.bio1, col="gray", legend=FALSE, main="toadflax Block")

points(block4, col="black", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))


#checkerboard

check <- get.checkerboard2(d[[1]], mean.temp.bio1, d[[2]], 25)

chk_blocks <- cbind(d[[1]], check$occ.grp)

block1 <- chk_blocks[which(chk_blocks[,3]==1),]

block2 <- chk_blocks[which(chk_blocks[,3]==2),]

block3 <- chk_blocks[which(chk_blocks[,3]==3),]

block4 <- chk_blocks[which(chk_blocks[,3]==4),]

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Toadflax Checkerboard")

points(block4, col="black", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))




#parsnip

bias_list <- get(paste(species.names[5], "_bias_list", sep="")) #get the species_bias_list object 

#prepare species occurrence and background data
d <- prepareData(species.file.paths[5], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

blocks <- get.block(d[[1]], d[[2]])

occ_blocks <- cbind(d[[1]], blocks$occ.grp)

block1 <- occ_blocks[which(occ_blocks[,3]==1),]

block2 <- occ_blocks[which(occ_blocks[,3]==2),]

block3 <- occ_blocks[which(occ_blocks[,3]==3),]

block4 <- occ_blocks[which(occ_blocks[,3]==4),]

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Parsnip Block")

points(block4, col="grey", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))


#checkerboard

check <- get.checkerboard2(d[[1]], mean.temp.bio1, d[[2]], 5)

chk_blocks <- cbind(d[[1]], check$occ.grp)

block1 <- chk_blocks[which(chk_blocks[,3]==1),]

block2 <- chk_blocks[which(chk_blocks[,3]==2),]

block3 <- chk_blocks[which(chk_blocks[,3]==3),]

block4 <- chk_blocks[which(chk_blocks[,3]==4),]

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Parsnip Checkerboard")

points(block4, col="grey", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))




#tansy

bias_list <- get(paste(species.names[6], "_bias_list", sep="")) #get the species_bias_list object 

#prepare species occurrence and background data
d <- prepareData(species.file.paths[6], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

#block

blocks <- get.block(d[[1]], d[[2]])

occ_blocks <- cbind(d[[1]], blocks$occ.grp)

block1 <- occ_blocks[which(occ_blocks[,3]==1),]

block2 <- occ_blocks[which(occ_blocks[,3]==2),]

block3 <- occ_blocks[which(occ_blocks[,3]==3),]

block4 <- occ_blocks[which(occ_blocks[,3]==4),]

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Tansy Block")

points(block4, col="grey", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))

#checkerboard

check <- get.checkerboard2(d[[1]], mean.temp.bio1, d[[2]], 5)

chk_blocks <- cbind(d[[1]], check$occ.grp)

block1 <- chk_blocks[which(chk_blocks[,3]==1),]

block2 <- chk_blocks[which(chk_blocks[,3]==2),]

block3 <- chk_blocks[which(chk_blocks[,3]==3),]

block4 <- chk_blocks[which(chk_blocks[,3]==4),]

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Tansy Checkerboard")

points(block4, col="grey", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))



#bittercress

bias_list <- get(paste(species.names[8], "_bias_list", sep="")) #get the species_bias_list object 

#prepare species occurrence and background data
d <- prepareData(species.file.paths[8], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

#block

blocks <- get.block(d[[1]], d[[2]])

occ_blocks <- cbind(d[[1]], blocks$occ.grp)

block1 <- occ_blocks[which(occ_blocks[,3]==1),]

block2 <- occ_blocks[which(occ_blocks[,3]==2),]

block3 <- occ_blocks[which(occ_blocks[,3]==3),]

block4 <- occ_blocks[which(occ_blocks[,3]==4),]

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Bittercress Block")

points(block4, col="grey", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))

#checkerboard

check <- get.checkerboard2(d[[1]], mean.temp.bio1, d[[2]], 5)

chk_blocks <- cbind(d[[1]], check$occ.grp)

block1 <- chk_blocks[which(chk_blocks[,3]==1),]

block2 <- chk_blocks[which(chk_blocks[,3]==2),]

block3 <- chk_blocks[which(chk_blocks[,3]==3),]

block4 <- chk_blocks[which(chk_blocks[,3]==4),]

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Tansy Checkerboard")

points(block4, col="grey", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))




#japanese hops

bias_list <- get(paste(species.names[7], "_bias_list", sep="")) #get the species_bias_list object 

#prepare species occurrence and background data
d <- prepareData(species.file.paths[7], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

#block

blocks <- get.block(d[[1]], d[[2]])

occ_blocks <- cbind(d[[1]], blocks$occ.grp)

block1 <- occ_blocks[which(occ_blocks[,3]==1),]

block2 <- occ_blocks[which(occ_blocks[,3]==2),]

block3 <- occ_blocks[which(occ_blocks[,3]==3),]

block4 <- occ_blocks[which(occ_blocks[,3]==4),]

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Japanese Hops Block")

points(block4, col="grey", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))

#checkerboard

check <- get.checkerboard2(d[[1]], mean.temp.bio1, d[[2]], 5)

chk_blocks <- cbind(d[[1]], check$occ.grp)

block1 <- chk_blocks[which(chk_blocks[,3]==1),]

block2 <- chk_blocks[which(chk_blocks[,3]==2),]

block3 <- chk_blocks[which(chk_blocks[,3]==3),]

block4 <- chk_blocks[which(chk_blocks[,3]==4),]

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Tansy Checkerboard")

points(block4, col="grey", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))


#oriental bittersweet

bias_list <- get(paste(species.names[9], "_bias_list", sep="")) #get the species_bias_list object 

#prepare species occurrence and background data
d <- prepareData(species.file.paths[9], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

#block

blocks <- get.block(d[[1]], d[[2]])

occ_blocks <- cbind(d[[1]], blocks$occ.grp)

block1 <- occ_blocks[which(occ_blocks[,3]==1),]

block2 <- occ_blocks[which(occ_blocks[,3]==2),]

block3 <- occ_blocks[which(occ_blocks[,3]==3),]

block4 <- occ_blocks[which(occ_blocks[,3]==4),]

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Oriental Bittersweet Block")

points(block4, col="grey", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))

#checkerboard

check <- get.checkerboard2(d[[1]], mean.temp.bio1, d[[2]], 5)

chk_blocks <- cbind(d[[1]], check$occ.grp)

block1 <- chk_blocks[which(chk_blocks[,3]==1),]

block2 <- chk_blocks[which(chk_blocks[,3]==2),]

block3 <- chk_blocks[which(chk_blocks[,3]==3),]

block4 <- chk_blocks[which(chk_blocks[,3]==4),]

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Tansy Checkerboard")

points(block4, col="grey", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))





















#Current Model Predictions

#####################

#Testing
setwd('c:/users/thomas/desktop/ENMEval Models')
folders <- list.files()


base <- "c:/users/thomas/desktop/ENMEval Models/"
fp <- paste(base, folders[1], sep="")
rda_files <- list.files(fp, full.names = TRUE)


for(j in 1:length(rda_files)){
  
  rda_file <- load(rda_files[j])
  
  cat(sprintf("\n------- Running Prediction %s -------\n", j))
  
  maxent.prediction <- rmaxent::project(best.model, bioclim.nb.3)
  
  #plot(maxent.prediction$prediction_cloglog)
  
  save(maxent.prediction, file = paste(tools::file_path_sans_ext(rda_files[j]), "_prediction.rda", sep=""))

  rm(rda_file, best.model, eval, thresholds, thresh.index, confusion, tss, kappa, boyce)
  
}


base <- "c:/users/thomas/desktop/ENMEval Models/"
fp <- paste(base, folders[2], sep="")
rda_files <- list.files(fp, full.names = TRUE)


for(j in 1:length(rda_files)){
  
  rda_file <- load(rda_files[j])
  
  cat(sprintf("\n------- Running Prediction %s -------\n", j))
  
  maxent.prediction <- rmaxent::project(best.model, bioclim.jh.3)
  
  #plot(maxent.prediction$prediction_cloglog)
  
  save(maxent.prediction, file = paste(tools::file_path_sans_ext(rda_files[j]), "_prediction.rda", sep=""))
  
  rm(rda_file, best.model, eval, thresholds, thresh.index, confusion, tss, kappa, boyce)
  
}



base <- "c:/users/thomas/desktop/ENMEval Models/"
fp <- paste(base, folders[3], sep="")
rda_files <- list.files(fp, full.names = TRUE)


for(j in 1:length(rda_files)){
  
  rda_file <- load(rda_files[j])
  
  cat(sprintf("\n------- Running Prediction %s -------\n", j))
  
  maxent.prediction <- rmaxent::project(best.model, bioclim.ob.3)
  
  #plot(maxent.prediction$prediction_cloglog)
  
  save(maxent.prediction, file = paste(tools::file_path_sans_ext(rda_files[j]), "_prediction.rda", sep=""))
  
  rm(rda_file, best.model, eval, thresholds, thresh.index, confusion, tss, kappa, boyce)
  
}



#####################



#Future Model Predictions

setwd('c:/users/thomas/desktop/ENMEval Models')
folders <- list.files()


base <- "c:/users/thomas/desktop/ENMEval Models/"
fp <- paste(base, folders[3], sep="")
rda_files <- list.files(fp, full.names = TRUE)

for(j in 1:length(rda_files)){

rda_file <- load(rda_files[j])

cat(sprintf("\n------- Running Future Predictions on Model %s -------\n", j))

start.time <- Sys.time()

#Predict Future Distribuitons
future.prediction.ccsm4_50 <- rmaxent::project(best.model, climate.stack.ccsm4_50)
save(future.prediction.ccsm4_50, file = paste(tools::file_path_sans_ext(rda_files[j]), "_ccsm4_50.rda", sep=""))

future.prediction.ccsm4_70 <- rmaxent::project(best.model, climate.stack.ccsm4_70)
save(future.prediction.ccsm4_70, file = paste(tools::file_path_sans_ext(rda_files[j]), "_ccsm4_70.rda", sep=""))

future.prediction.gfdl_50 <- rmaxent::project(best.model, climate.stack.gfdl_50)
save(future.prediction.gfdl_50, file = paste(tools::file_path_sans_ext(rda_files[j]), "_gfdl_50.rda", sep=""))

future.prediction.gfdl_70 <- rmaxent::project(best.model, climate.stack.gfdl_70)
save(future.prediction.gfdl_70, file = paste(tools::file_path_sans_ext(rda_files[j]), "_gfdl_70.rda", sep=""))

future.prediction.ipsl_50 <- rmaxent::project(best.model, climate.stack.ipsl_50)
save(future.prediction.ipsl_50, file = paste(tools::file_path_sans_ext(rda_files[j]), "_ipsl_50.rda", sep=""))

future.prediction.ipsl_70 <- rmaxent::project(best.model, climate.stack.ipsl_70)
save(future.prediction.ipsl_70, file = paste(tools::file_path_sans_ext(rda_files[j]), "_ipsl_70.rda", sep=""))

future.prediction.mri_50 <- rmaxent::project(best.model, climate.stack.mri_50)
save(future.prediction.mri_50, file = paste(tools::file_path_sans_ext(rda_files[j]), "_mri_50.rda", sep=""))

future.prediction.mri_70 <- rmaxent::project(best.model, climate.stack.mri_70)
save(future.prediction.mri_70, file = paste(tools::file_path_sans_ext(rda_files[j]), "_mri_70.rda", sep=""))

future.prediction.miroc_50 <- rmaxent::project(best.model, climate.stack.miroc_50)
save(future.prediction.miroc_50, file = paste(tools::file_path_sans_ext(rda_files[j]), "_miroc_50.rda", sep=""))

future.prediction.miroc_70 <- rmaxent::project(best.model, climate.stack.miroc_70)
save(future.prediction.miroc_70, file = paste(tools::file_path_sans_ext(rda_files[j]), "_miroc_70.rda", sep=""))

end.time <- Sys.time()

print(end.time - start.time)

gc()

rm(rda_file, best.model, eval, thresholds, thresh.index, confusion, tss, kappa, boyce)


}


#####



#Plots for examining downsampling effects


us.coord.extent <- extent(-135, -50, 8.8, 60)
mw.coord.extent <- extent(-100, -80, 35, 50)
east.coord.extent <- extent(-74, -60, 8.8, 45)

bioclim.crop <- crop(bioclim.ct.3, mw.coord.extent)

file <- sw.file.path
species.data <- read.csv(file, header=TRUE, sep=",")
species.coords <- species.data[,2:3]
coordinates(species.coords) <- ~ lon + lat
projection(species.coords) <- CRS('+proj=longlat +datum=WGS84')
plot(species.coords) #162 species occurrence records

species.coords.east <- crop(species.coords, east.coord.extent)
plot(species.coords.east)
points.0km = length(species.coords.east)

sample.grid.mask <- mean.temp.bio1
res(sample.grid.mask) <- 0.001 #create a raster with grid cells of 0.008 degrees
sample.grid.mask <- extend(sample.grid.mask, extent(sample.grid.mask)+0.001)
subsamp.species.coords.1km <- gridSample(species.coords.east, sample.grid.mask, n = 1)
plot(subsamp.species.coords.1km)
points.1km = length(subsamp.species.coords.1km)

sample.grid.mask <- mean.temp.bio1
res(sample.grid.mask) <- 0.01 #create a raster with grid cells of 0.008 degrees
sample.grid.mask <- extend(sample.grid.mask, extent(sample.grid.mask)+0.01)
subsamp.species.coords.10km <- gridSample(species.coords.east, sample.grid.mask, n = 1)
plot(subsamp.species.coords.10km)
points.10km = length(subsamp.species.coords.10km)



#Plots for examining GDK effects





#Plots for testing model questions


#toadflax
setwd('c:/users/thomas/desktop/ENMEval Models')
folders <- list.files()

base <- "c:/users/thomas/desktop/ENMEval Models/"
fp <- paste(base, folders[9], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#max train AUC
rda_file <- load(rda_files[7])

maxent.prediction_max_auc <- rmaxent::project(best.model, bioclim.tf.3)

plot(maxent.prediction_max_auc$prediction_cloglog, main = "Toadflax Max AUC Model")

#max auc aicc 0
rda_file <- load(rda_files[13])

maxent.prediction_delta0_max_auc <- rmaxent::project(best.model, bioclim.tf.3)

plot(maxent.prediction_delta0_max_auc$prediction_cloglog, main = "Toadflax Max AUC Model with DeltaAIC = 0")

#min var model
rda_file <- load(rda_files[14])

maxent.prediction_min_var <- rmaxent::project(best.model, bioclim.tf.3)

plot(maxent.prediction_min_var$prediction_cloglog, main = "Toadflax Min Var with DeltaAIC = 0")


#plot 3 models for comparsion
par(mfrow=c(3,1))

plot(maxent.prediction_max_auc$prediction_cloglog, main = "Toadflax Max AUC Model")

plot(maxent.prediction_delta0_max_auc$prediction_cloglog, main = "Toadflax Max AUC Model with DeltaAIC = 0")

plot(maxent.prediction_min_var$prediction_cloglog, main = "Toadflax Min Var with DeltaAIC = 0")










#knapweed
setwd('c:/users/thomas/desktop/ENMEval Models')
folders <- list.files()

base <- "c:/users/thomas/desktop/ENMEval Models/"
fp <- paste(base, folders[4], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#max train AUC
rda_file <- load(rda_files[13])

maxent.prediction_max_auc <- rmaxent::project(best.model, bioclim.bk.3)

#max auc aicc 0
rda_file <- load(rda_files[13])

maxent.prediction_delta0_max_auc <- rmaxent::project(best.model, bioclim.bk.3)

#min var model
rda_file <- load(rda_files[14])

maxent.prediction_min_var <- rmaxent::project(best.model, bioclim.bk.3)

#plot 3 models for comparsion
par(mfrow=c(3,1))

plot(maxent.prediction_max_auc$prediction_cloglog, main = "Knapweed Max AUC Model")

plot(maxent.prediction_delta0_max_auc$prediction_cloglog, main = "Knapweed Max AUC Model with DeltaAIC = 0")

plot(maxent.prediction_min_var$prediction_cloglog, main = "Knapweed Min Var with DeltaAIC = 0")






#parsnip
setwd('c:/users/thomas/desktop/ENMEval Models')
folders <- list.files()

base <- "c:/users/thomas/desktop/ENMEval Models/"
fp <- paste(base, folders[5], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#max train AUC
rda_file <- load(rda_files[22])

maxent.prediction_max_auc <- rmaxent::project(best.model, bioclim.parsnip.3)

#max auc aicc 0
rda_file <- load(rda_files[23])

maxent.prediction_delta0_max_auc <- rmaxent::project(best.model, bioclim.parsnip.3)

#min var model
rda_file <- load(rda_files[15])

maxent.prediction_min_var <- rmaxent::project(best.model, bioclim.parsnip.3)

#plot 3 models for comparsion
par(mfrow=c(3,1))

plot(maxent.prediction_max_auc$prediction_cloglog, main = "Parsnip Max AUC Model")

plot(maxent.prediction_delta0_max_auc$prediction_cloglog, main = "Parsnip Max AUC Model with DeltaAIC = 0")

plot(maxent.prediction_min_var$prediction_cloglog, main = "Parsnip Min Var with DeltaAIC = 0")






#swallowwort
setwd('c:/users/thomas/desktop/ENMEval Models')
folders <- list.files()

base <- "c:/users/thomas/desktop/ENMEval Models/"
fp <- paste(base, folders[6], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#max train AUC
rda_file <- load(rda_files[1])

maxent.prediction_max_auc <- rmaxent::project(best.model, bioclim.sw.3)

#max auc aicc 0
rda_file <- load(rda_files[5])

maxent.prediction_delta0_max_auc <- rmaxent::project(best.model, bioclim.sw.3)

#min var model
rda_file <- load(rda_files[6])

maxent.prediction_min_var <- rmaxent::project(best.model, bioclim.sw.3)

#plot 3 models for comparsion
par(mfrow=c(3,1))

plot(maxent.prediction_max_auc$prediction_cloglog, main = "Swallowwort Max AUC Model")

plot(maxent.prediction_delta0_max_auc$prediction_cloglog, main = "Swallowwort Max AUC Model with DeltaAIC = 0")

plot(maxent.prediction_min_var$prediction_cloglog, main = "Swallowwort Min Var with DeltaAIC = 0")





#tansy
setwd('c:/users/thomas/desktop/ENMEval Models')
folders <- list.files()

base <- "c:/users/thomas/desktop/ENMEval Models/"
fp <- paste(base, folders[7], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#max train AUC
rda_file <- load(rda_files[3])

maxent.prediction_max_auc <- rmaxent::project(best.model, bioclim.tansy.3)

#max auc aicc 0
rda_file <- load(rda_files[3])

maxent.prediction_delta0_max_auc <- rmaxent::project(best.model, bioclim.tansy.3)

#min var model
rda_file <- load(rda_files[4])

maxent.prediction_min_var <- rmaxent::project(best.model, bioclim.tansy.3)

#plot 3 models for comparsion
par(mfrow=c(3,1))

plot(maxent.prediction_max_auc$prediction_cloglog, main = "Tansy Max AUC Model")

plot(maxent.prediction_delta0_max_auc$prediction_cloglog, main = "Tansy Max AUC Model with DeltaAIC = 0")

plot(maxent.prediction_min_var$prediction_cloglog, main = "Tansy Min Var with DeltaAIC = 0")





#teasel
setwd('c:/users/thomas/desktop/ENMEval Models')
folders <- list.files()

base <- "c:/users/thomas/desktop/ENMEval Models/"
fp <- paste(base, folders[8], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#max train AUC
rda_file <- load(rda_files[1])

maxent.prediction_max_auc <- rmaxent::project(best.model, bioclim.ct.3)

#max auc aicc 0
rda_file <- load(rda_files[2])

maxent.prediction_delta0_max_auc <- rmaxent::project(best.model, bioclim.ct.3)

#min var model
rda_file <- load(rda_files[21])

maxent.prediction_min_var <- rmaxent::project(best.model, bioclim.ct.3)

#plot 3 models for comparsion
par(mfrow=c(3,1))

plot(maxent.prediction_max_auc$prediction_cloglog, main = "Teasel Max AUC Model")

plot(maxent.prediction_delta0_max_auc$prediction_cloglog, main = "Teasel Max AUC Model with DeltaAIC = 0")

plot(maxent.prediction_min_var$prediction_cloglog, main = "Teasel Min Var with DeltaAIC = 0")





















####################
# Function Testing #
####################


#file_base <- sprintf("F:/ENMEval Models Reduced Parameter Set 5-1-2019/%s/", species.names[1])

source('c:/users/thomas/desktop/enmeval_functions.R')

d <- prepareData(species.file.paths[2], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

o <- runENMModel(d[[1]], env.predictors[[2]], d[[2]], species_name = species.names[2], gdk_value = gdk.values[1], downsample_value = downsample.units[1])






eval_model <- readRDS("F:/ENMEval Models Reduced Parameter Set 5-1-2019/teasel/species_teasel_gdk_0_ds_1km_method_block_full_eval_model.rds")

print(eval_model@results)
#select "best" model, which has lowest delta AIC (likely 0)
#model.index <- which(eval_model@results$delta.AICc==0)

#Error in model index as an incrementor for the models

model.index <- NROW(eval_model@results)

print(model.index)

for(i in 1:length(model.index)){ #iterate through all models where delta AIC=0
  
  best.model <- eval_model@models[[model.index[i]]]
  #best.model <- eval_model@models[[model.index[1]]]#in the case of multiple models with AIC=0, select all models with AIC=0
  
  print(best.model)
  
  #calculate all model statistics
  library(dismo)
  eval <- dismo::evaluate(p=subsampled_species_coordinates, a=background_coordinates, model=eval_model@models[[1]], x=environment)
  print(eval)
  thresholds <- threshold(eval)
  thresh.index <- which((eval@TPR+eval@TNR) == max(eval@TPR+eval@TNR))
  kappa <- eval@kappa[thresh.index[1]]
  #wrong kappa
  confusion <- eval@confusion[thresh.index,]
  print(confusion)
  tss<- max(eval@TPR+eval@TNR)-1
  #wrong boyce
  boyce <- contBoyce2x(pres=eval@presence, bg=eval@absence, numClasses = 10, upweightTails = TRUE, na.rm=TRUE, autoWindow = TRUE, graph = TRUE, method = 'pearson')
  
  #initialize model name for output
  #model_name <- paste(file_base, sprintf("species_%s_gdk_%s_ds_%s_method_%s_model_%s.rda", species.names[1], gdk.values[1], downsample.values[1], "block method", model.index[i]), sep="") #model name goes here for file name output
  
  #print(model_name)
  
  #save model output
  #save(best.model, eval, thresholds, thresh.index, confusion, tss, kappa, boyce, file = model_name, compress = TRUE)
 # rm(best.model, eval, thresholds, thresh.index, confusion, tss, kappa, boyce)
}
#}

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  


  
  #################################################################################
  
#Large results, master figure to plot species Block, Checkerboard "Best" Models
  
####################################################################################
  
setwd("F:/ENMEval Models Reduced Parameter Set 5-1-2019")


#BITTERCRESS


bias_list <- get(paste(species.names[8], "_bias_list", sep="")) #get the species_bias_list object 

d <- prepareData(species.file.paths[8], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[])


folders <- list.files()

base <- "F:/ENMEval Models Reduced Parameter Set 5-1-2019/"
fp <- paste(base, folders[1], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#block

#max test AUC
rda_file <- load(rda_files[12])

maxent.prediction_max_auc_b <- rmaxent::project(best.model, bioclim.nb.3)

#min aicc
rda_file <- load(rda_files[5])

maxent.prediction_min_aicc_b <- rmaxent::project(best.model, bioclim.nb.3)

#min var between validation groups
rda_file <- load(rda_files[2])

maxent.prediction_min_var_b <- rmaxent::project(best.model, bioclim.nb.3)

#max sensitivity
rda_file <- load(rda_files[45])

maxent.prediction_max_sens_b <- rmaxent::project(best.model, bioclim.nb.3)

#checkerboard

#max test AUC
rda_file <- load(rda_files[17])

maxent.prediction_max_auc_c <- rmaxent::project(best.model, bioclim.nb.3)

#min aicc
rda_file <- load(rda_files[10])

maxent.prediction_min_aicc_c <- rmaxent::project(best.model, bioclim.nb.3)

#min var between validation groups
rda_file <- load(rda_files[7])

maxent.prediction_min_var_c <- rmaxent::project(best.model, bioclim.nb.3)

#max sensitivity
rda_file <- load(rda_files[48])

maxent.prediction_max_sens_c <- rmaxent::project(best.model, bioclim.nb.3)



#block row
par(mfrow=c(2,4))
plot(maxent.prediction_max_auc_b$prediction_cloglog, main = "Bittercress Max AUC Model Block")

plot(maxent.prediction_min_aicc_b$prediction_cloglog, main = "Bittercress Min AICc Model Block")


plot(maxent.prediction_min_var_b$prediction_cloglog, main = "Bittercress Min Var Model Block")



plot(maxent.prediction_max_sens_b$prediction_cloglog, main = "Bittercress Max Sens Model Block")
#checkerboard row

plot(maxent.prediction_max_auc_c$prediction_cloglog, main = "Bittercress Max AUC Model Checkerboard")



plot(maxent.prediction_min_aicc_c$prediction_cloglog, main = "Bittercress Min AICc Model Checkerboard")



plot(maxent.prediction_min_var_c$prediction_cloglog, main = "Bittercress Min Var Model Checkerboard")


plot(maxent.prediction_max_sens_c$prediction_cloglog, main = "Bittercress Max Sens Model Checkerboard")







#BITTERSWEET


bias_list <- get(paste(species.names[9], "_bias_list", sep="")) #get the species_bias_list object 

d <- prepareData(species.file.paths[9], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[])


folders <- list.files()

base <- "F:/ENMEval Models Reduced Parameter Set 5-1-2019/"
fp <- paste(base, folders[2], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#block

#max test AUC
rda_file <- load(rda_files[15])

maxent.prediction_max_auc_b <- rmaxent::project(best.model, bioclim.ob.3)

#min aicc
rda_file <- load(rda_files[2])

maxent.prediction_min_aicc_b <- rmaxent::project(best.model, bioclim.ob.3)

#min var between validation groups
rda_file <- load(rda_files[35])

maxent.prediction_min_var_b <- rmaxent::project(best.model, bioclim.ob.3)

#max sensitivity
rda_file <- load(rda_files[14])

maxent.prediction_max_sens_b <- rmaxent::project(best.model, bioclim.ob.3)

#checkerboard

#max test AUC
rda_file <- load(rda_files[17])

maxent.prediction_max_auc_c <- rmaxent::project(best.model, bioclim.ob.3)

#min aicc
rda_file <- load(rda_files[7])

maxent.prediction_min_aicc_c <- rmaxent::project(best.model, bioclim.ob.3)

#min var between validation groups
rda_file <- load(rda_files[17])

maxent.prediction_min_var_c <- rmaxent::project(best.model, bioclim.ob.3)

#max sensitivity
rda_file <- load(rda_files[19])

maxent.prediction_max_sens_c <- rmaxent::project(best.model, bioclim.ob.3)


#plot 8 models for comparsion
par(mfrow=c(2,4))
#block row

plot(maxent.prediction_max_auc_b$prediction_cloglog, main = "bittersweet Max AUC Model Block")

plot(maxent.prediction_min_aicc_b$prediction_cloglog, main = "bittersweet Min AICc Model Block")

plot(maxent.prediction_min_var_b$prediction_cloglog, main = "bittersweet Min Var Model Block")


plot(maxent.prediction_max_sens_b$prediction_cloglog, main = "bittersweet Max Sens Model Block")
#checkerboard row
plot(maxent.prediction_max_auc_c$prediction_cloglog, main = "bittersweet Max AUC Model Checkerboard")


plot(maxent.prediction_min_aicc_c$prediction_cloglog, main = "bittersweet Min AICc Model Checkerboard")


plot(maxent.prediction_min_var_c$prediction_cloglog, main = "bittersweet Min Var Model Checkerboard")


plot(maxent.prediction_max_sens_c$prediction_cloglog, main = "bittersweet Max Sens Model Checkerboard")





#Knapweed

bias_list <- get(paste(species.names[3], "_bias_list", sep="")) #get the species_bias_list object 

d <- prepareData(species.file.paths[3], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[])


folders <- list.files()

base <- "F:/ENMEval Models Reduced Parameter Set 5-1-2019/"
fp <- paste(base, folders[4], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#block

#max test AUC
rda_file <- load(rda_files[2])

maxent.prediction_max_auc_b <- rmaxent::project(best.model, bioclim.bk.3)

#min aicc
rda_file <- load(rda_files[3])

maxent.prediction_min_aicc_b <- rmaxent::project(best.model, bioclim.bk.3)

#min var between validation groups
rda_file <- load(rda_files[2])

maxent.prediction_min_var_b <- rmaxent::project(best.model, bioclim.bk.3)

#max sensitivity
rda_file <- load(rda_files[4])

maxent.prediction_max_sens_b <- rmaxent::project(best.model, bioclim.bk.3)

#checkerboard

#max test AUC
rda_file <- load(rda_files[17])

maxent.prediction_max_auc_c <- rmaxent::project(best.model, bioclim.bk.3)

#min aicc
rda_file <- load(rda_files[8])

maxent.prediction_min_aicc_c <- rmaxent::project(best.model, bioclim.bk.3)

#min var between validation groups
rda_file <- load(rda_files[30])

maxent.prediction_min_var_c <- rmaxent::project(best.model, bioclim.bk.3)

#max sensitivity
rda_file <- load(rda_files[10])

maxent.prediction_max_sens_c <- rmaxent::project(best.model, bioclim.bk.3)


#plot 8 models for comparsion
par(mfrow=c(2,4))
#block row
plot(maxent.prediction_max_auc_b$prediction_cloglog, main = "Knapweed Max AUC Model Block")

plot(maxent.prediction_min_aicc_b$prediction_cloglog, main = "Knapweed Min AICc Model Block")

plot(maxent.prediction_min_var_b$prediction_cloglog, main = "Knapweed Min Var Model Block")

plot(maxent.prediction_max_sens_b$prediction_cloglog, main = "Knapweed Max Sens Model Block")
#checkerboard row
plot(maxent.prediction_max_auc_c$prediction_cloglog, main = "Knapweed Max AUC Model Checkerboard")

plot(maxent.prediction_min_aicc_c$prediction_cloglog, main = "Knapweed Min AICc Model Checkerboard")

plot(maxent.prediction_min_var_c$prediction_cloglog, main = "Knapweed Min Var Model Checkerboard")

plot(maxent.prediction_max_sens_c$prediction_cloglog, main = "Knapweed Max Sens Model Checkerboard")




#Japanese Hops

bias_list <- get(paste(species.names[7], "_bias_list", sep="")) #get the species_bias_list object 

d <- prepareData(species.file.paths[7], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[])


folders <- list.files()

base <- "F:/ENMEval Models Reduced Parameter Set 5-1-2019/"
fp <- paste(base, folders[3], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#block

#max test AUC

rda_file <- load(rda_files[5])

maxent.prediction_max_auc_b <- rmaxent::project(best.model, bioclim.jh.3)

#min aicc
rda_file <- load(rda_files[2])

maxent.prediction_min_aicc_b <- rmaxent::project(best.model, bioclim.jh.3)

#min var between validation groups
rda_file <- load(rda_files[2])

maxent.prediction_min_var_b <- rmaxent::project(best.model, bioclim.jh.3)

#max sensitivity
rda_file <- load(rda_files[5])

maxent.prediction_max_sens_b <- rmaxent::project(best.model, bioclim.jh.3)

#checkerboard

#max test AUC
rda_file <- load(rda_files[18])

maxent.prediction_max_auc_c <- rmaxent::project(best.model, bioclim.jh.3)

#min aicc
rda_file <- load(rda_files[7])

maxent.prediction_min_aicc_c <- rmaxent::project(best.model, bioclim.jh.3)

#min var between validation groups
rda_file <- load(rda_files[7])

maxent.prediction_min_var_c <- rmaxent::project(best.model, bioclim.jh.3)

#max sensitivity
rda_file <- load(rda_files[10])

maxent.prediction_max_sens_c <- rmaxent::project(best.model, bioclim.jh.3)


#plot 8 models for comparsion
par(mfrow=c(2,4))
#block row
plot(maxent.prediction_max_auc_b$prediction_cloglog, main = "Japanese Hops Max AUC Model Block")

plot(maxent.prediction_min_aicc_b$prediction_cloglog, main = "Japanese Hops Min AICc Model Block")

plot(maxent.prediction_min_var_b$prediction_cloglog, main = "Japanese Hops Min Var Model Block")

plot(maxent.prediction_max_sens_b$prediction_cloglog, main = "Japanese Hops Max Sens Model Block")
#checkerboard row
plot(maxent.prediction_max_auc_c$prediction_cloglog, main = "Japanese Hops Max AUC Model Checkerboard")

plot(maxent.prediction_min_aicc_c$prediction_cloglog, main = "Japanese Hops Min AICc Model Checkerboard")

plot(maxent.prediction_min_var_c$prediction_cloglog, main = "Japanese Hops Min Var Model Checkerboard")

plot(maxent.prediction_max_sens_c$prediction_cloglog, main = "Japanese Hops Max Sens Model Checkerboard")




#parsnip

bias_list <- get(paste(species.names[5], "_bias_list", sep="")) #get the species_bias_list object 

d <- prepareData(species.file.paths[5], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[])


folders <- list.files()

base <- "F:/ENMEval Models Reduced Parameter Set 5-1-2019/"
fp <- paste(base, folders[5], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#block

#max test AUC

rda_file <- load(rda_files[13])

maxent.prediction_max_auc_b <- rmaxent::project(best.model, bioclim.parsnip.3)

#min aicc
rda_file <- load(rda_files[2])

maxent.prediction_min_aicc_b <- rmaxent::project(best.model, bioclim.parsnip.3)

#min var between validation groups
rda_file <- load(rda_files[13])

maxent.prediction_min_var_b <- rmaxent::project(best.model, bioclim.parsnip.3)

#max sensitivity
rda_file <- load(rda_files[5])

maxent.prediction_max_sens_b <- rmaxent::project(best.model, bioclim.parsnip.3)

#checkerboard

#max test AUC
rda_file <- load(rda_files[12])

maxent.prediction_max_auc_c <- rmaxent::project(best.model, bioclim.parsnip.3)

#min aicc
rda_file <- load(rda_files[7])

maxent.prediction_min_aicc_c <- rmaxent::project(best.model, bioclim.parsnip.3)

#min var between validation groups
rda_file <- load(rda_files[58])

maxent.prediction_min_var_c <- rmaxent::project(best.model, bioclim.parsnip.3)

#max sensitivity
rda_file <- load(rda_files[17])

maxent.prediction_max_sens_c <- rmaxent::project(best.model, bioclim.parsnip.3)


#plot 8 models for comparsion
par(mfrow=c(2,4))
#block row
plot(maxent.prediction_max_auc_b$prediction_cloglog, main = "Parsnip Max AUC Model Block")

plot(maxent.prediction_min_aicc_b$prediction_cloglog, main = "Parsnip Min AICc Model Block")

plot(maxent.prediction_min_var_b$prediction_cloglog, main = "Parsnip Min Var Model Block")

plot(maxent.prediction_max_sens_b$prediction_cloglog, main = "Parsnip Max Sens Model Block")
#checkerboard row
plot(maxent.prediction_max_auc_c$prediction_cloglog, main = "Parsnip Max AUC Model Checkerboard")

plot(maxent.prediction_min_aicc_c$prediction_cloglog, main = "Parsnip Min AICc Model Checkerboard")

plot(maxent.prediction_min_var_c$prediction_cloglog, main = "Parsnip Min Var Model Checkerboard")

plot(maxent.prediction_max_sens_c$prediction_cloglog, main = "Parsnip Max Sens Model Checkerboard")





#swallowwort

bias_list <- get(paste(species.names[2], "_bias_list", sep="")) #get the species_bias_list object 

d <- prepareData(species.file.paths[2], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[])


folders <- list.files()

base <- "F:/ENMEval Models Reduced Parameter Set 5-1-2019/"
fp <- paste(base, folders[6], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#block

#max test AUC

rda_file <- load(rda_files[14])

maxent.prediction_max_auc_b <- rmaxent::project(best.model, bioclim.sw.3)

#min aicc
rda_file <- load(rda_files[3])

maxent.prediction_min_aicc_b <- rmaxent::project(best.model, bioclim.sw.3)

#min var between validation groups
rda_file <- load(rda_files[2])

maxent.prediction_min_var_b <- rmaxent::project(best.model, bioclim.sw.3)

#max sensitivity
rda_file <- load(rda_files[35])

maxent.prediction_max_sens_b <- rmaxent::project(best.model, bioclim.sw.3)

#checkerboard

#max test AUC
rda_file <- load(rda_files[17])

maxent.prediction_max_auc_c <- rmaxent::project(best.model, bioclim.sw.3)

#min aicc
rda_file <- load(rda_files[8])

maxent.prediction_min_aicc_c <- rmaxent::project(best.model, bioclim.sw.3)

#min var between validation groups
rda_file <- load(rda_files[17])

maxent.prediction_min_var_c <- rmaxent::project(best.model, bioclim.sw.3)

#max sensitivity
rda_file <- load(rda_files[40])

maxent.prediction_max_sens_c <- rmaxent::project(best.model, bioclim.sw.3)


#plot 8 models for comparsion
par(mfrow=c(2,4))
#block row
plot(maxent.prediction_max_auc_b$prediction_cloglog, main = "Swallowwort Max AUC Model Block")

plot(maxent.prediction_min_aicc_b$prediction_cloglog, main = "Swallowwort Min AICc Model Block")

plot(maxent.prediction_min_var_b$prediction_cloglog, main = "Swallowwort Min Var Model Block")

plot(maxent.prediction_max_sens_b$prediction_cloglog, main = "Swallowwort Max Sens Model Block")
#checkerboard row
plot(maxent.prediction_max_auc_c$prediction_cloglog, main = "Swallowwort Max AUC Model Checkerboard")

plot(maxent.prediction_min_aicc_c$prediction_cloglog, main = "Swallowwort Min AICc Model Checkerboard")

plot(maxent.prediction_min_var_c$prediction_cloglog, main = "Swallowwort Min Var Model Checkerboard")

plot(maxent.prediction_max_sens_c$prediction_cloglog, main = "Swallowwort Max Sens Model Checkerboard")





#tansy

bias_list <- get(paste(species.names[6], "_bias_list", sep="")) #get the species_bias_list object 

d <- prepareData(species.file.paths[6], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[])


folders <- list.files()

base <- "F:/ENMEval Models Reduced Parameter Set 5-1-2019/"
fp <- paste(base, folders[7], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#block

#max test AUC

rda_file <- load(rda_files[12])

maxent.prediction_max_auc_b <- rmaxent::project(best.model, bioclim.tansy.3)

#min aicc
rda_file <- load(rda_files[2])

maxent.prediction_min_aicc_b <- rmaxent::project(best.model, bioclim.tansy.3)

#min var between validation groups
rda_file <- load(rda_files[42])

maxent.prediction_min_var_b <- rmaxent::project(best.model, bioclim.tansy.3)

#max sensitivity
rda_file <- load(rda_files[53])

maxent.prediction_max_sens_b <- rmaxent::project(best.model, bioclim.tansy.3)

#checkerboard

#max test AUC
rda_file <- load(rda_files[17])

maxent.prediction_max_auc_c <- rmaxent::project(best.model, bioclim.tansy.3)

#min aicc
rda_file <- load(rda_files[7])

maxent.prediction_min_aicc_c <- rmaxent::project(best.model, bioclim.tansy.3)

#min var between validation groups
rda_file <- load(rda_files[17])

maxent.prediction_min_var_c <- rmaxent::project(best.model, bioclim.tansy.3)

#max sensitivity
#change this
rda_file <- load(rda_files[58])

maxent.prediction_max_sens_c <- rmaxent::project(best.model, bioclim.tansy.3)


#plot 8 models for comparsion
par(mfrow=c(2,4))
#block row
plot(maxent.prediction_max_auc_b$prediction_cloglog, main = "Common Tansy Max AUC Model Block")

plot(maxent.prediction_min_aicc_b$prediction_cloglog, main = "Common Tansy Min AICc Model Block")

plot(maxent.prediction_min_var_b$prediction_cloglog, main = "Common Tansy Min Var Model Block")

plot(maxent.prediction_max_sens_b$prediction_cloglog, main = "Common Tansy Max Sens Model Block")
#checkerboard row
plot(maxent.prediction_max_auc_c$prediction_cloglog, main = "Common Tansy Max AUC Model Checkerboard")

plot(maxent.prediction_min_aicc_c$prediction_cloglog, main = "Common Tansy Min AICc Model Checkerboard")

plot(maxent.prediction_min_var_c$prediction_cloglog, main = "Common Tansy Min Var Model Checkerboard")

plot(maxent.prediction_max_sens_c$prediction_cloglog, main = "Common Tansy Max Sens Model Checkerboard")







#teasel

bias_list <- get(paste(species.names[1], "_bias_list", sep="")) #get the species_bias_list object 

d <- prepareData(species.file.paths[1], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[])


folders <- list.files()

base <- "F:/ENMEval Models Reduced Parameter Set 5-1-2019/"
fp <- paste(base, folders[8], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#block

#max test AUC

rda_file <- load(rda_files[12])

maxent.prediction_max_auc_b <- rmaxent::project(best.model, bioclim.ct.3)

#min aicc
rda_file <- load(rda_files[3])

maxent.prediction_min_aicc_b <- rmaxent::project(best.model, bioclim.ct.3)

#min var between validation groups
rda_file <- load(rda_files[2])

maxent.prediction_min_var_b <- rmaxent::project(best.model, bioclim.ct.3)

#max sensitivity
rda_file <- load(rda_files[32])

maxent.prediction_max_sens_b <- rmaxent::project(best.model, bioclim.ct.3)

#checkerboard

#max test AUC
rda_file <- load(rda_files[17])

maxent.prediction_max_auc_c <- rmaxent::project(best.model, bioclim.ct.3)

#min aicc
rda_file <- load(rda_files[8])

maxent.prediction_min_aicc_c <- rmaxent::project(best.model, bioclim.ct.3)

#min var between validation groups
rda_file <- load(rda_files[27])

maxent.prediction_min_var_c <- rmaxent::project(best.model, bioclim.ct.3)

#max sensitivity
rda_file <- load(rda_files[37])

maxent.prediction_max_sens_c <- rmaxent::project(best.model, bioclim.ct.3)


#plot 8 models for comparsion
par(mfrow=c(2,4))
#block row
plot(maxent.prediction_max_auc_b$prediction_cloglog, main = "Common Teasel Max AUC Model Block")

plot(maxent.prediction_min_aicc_b$prediction_cloglog, main = "Common Teasel Min AICc Model Block")

plot(maxent.prediction_min_var_b$prediction_cloglog, main = "Common Teasel Min Var Model Block")

plot(maxent.prediction_max_sens_b$prediction_cloglog, main = "Common Teasel Max Sens Model Block")
#checkerboard row
plot(maxent.prediction_max_auc_c$prediction_cloglog, main = "Common Teasel Max AUC Model Checkerboard")

plot(maxent.prediction_min_aicc_c$prediction_cloglog, main = "Common Teasel Min AICc Model Checkerboard")

plot(maxent.prediction_min_var_c$prediction_cloglog, main = "Common Teasel Min Var Model Checkerboard")

plot(maxent.prediction_max_sens_c$prediction_cloglog, main = "Common Teasel Max Sens Model Checkerboard")




#toadflax

bias_list <- get(paste(species.names[4], "_bias_list", sep="")) #get the species_bias_list object 

d <- prepareData(species.file.paths[4], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[])


folders <- list.files()

base <- "F:/ENMEval Models Reduced Parameter Set 5-1-2019/"
fp <- paste(base, folders[9], sep="")
rda_files <- list.files(fp, full.names = TRUE)

#block

#max test AUC

rda_file <- load(rda_files[12])

maxent.prediction_max_auc_b <- rmaxent::project(best.model, bioclim.tf.3)

#min aicc
rda_file <- load(rda_files[2])

maxent.prediction_min_aicc_b <- rmaxent::project(best.model, bioclim.tf.3)

#min var between validation groups
rda_file <- load(rda_files[53])

maxent.prediction_min_var_b <- rmaxent::project(best.model, bioclim.tf.3)

#max sensitivity
rda_file <- load(rda_files[13])

maxent.prediction_max_sens_b <- rmaxent::project(best.model, bioclim.tf.3)

#checkerboard

#max test AUC
rda_file <- load(rda_files[17])

maxent.prediction_max_auc_c <- rmaxent::project(best.model, bioclim.tf.3)

#min aicc
rda_file <- load(rda_files[7])

maxent.prediction_min_aicc_c <- rmaxent::project(best.model, bioclim.tf.3)

#min var between validation groups
rda_file <- load(rda_files[40])

maxent.prediction_min_var_c <- rmaxent::project(best.model, bioclim.tf.3)

#max sensitivity
rda_file <- load(rda_files[18])

maxent.prediction_max_sens_c <- rmaxent::project(best.model, bioclim.tf.3)


#plot 8 models for comparsion
par(mfrow=c(2,4))
#block row
plot(maxent.prediction_max_auc_b$prediction_cloglog, main = "Toadflax Max AUC Model Block")

plot(maxent.prediction_min_aicc_b$prediction_cloglog, main = "Toadflax Min AICc Model Block")

plot(maxent.prediction_min_var_b$prediction_cloglog, main = "Toadflax Min Var Model Block")

plot(maxent.prediction_max_sens_b$prediction_cloglog, main = "Toadflax Max Sens Model Block")
#checkerboard row
plot(maxent.prediction_max_auc_c$prediction_cloglog, main = "Toadflax Max AUC Model Checkerboard")

plot(maxent.prediction_min_aicc_c$prediction_cloglog, main = "Toadflax Min AICc Model Checkerboard")

plot(maxent.prediction_min_var_c$prediction_cloglog, main = "Toadflax Min Var Model Checkerboard")

plot(maxent.prediction_max_sens_c$prediction_cloglog, main = "Toadflax Max Sens Model Checkerboard")














####################################

#Plotting Results May 15 2019

####################################


######################################################

setwd('F:/ENMEval Models')

folders <-list.files()

#######################
# BITTERCRESS RESULTS #
#######################

base <- "F:/ENMEval Models/"
fp <- paste(base, folders[1], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#gdk 0 beta 1
rda_file <- load(rda_files[1])
bittercress_gdk0b1 <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 2
rda_file <- load(rda_files[2])
bittercress_gdk0b2 <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 3
rda_file <- load(rda_files[3])
bittercress_gdk0b3 <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 1
rda_file <- load(rda_files[17])
bittercress_gdk1b1 <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 2
rda_file <- load(rda_files[18])
bittercress_gdk1b2 <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 3
rda_file <- load(rda_files[19])
bittercress_gdk1b3 <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 1
rda_file <- load(rda_files[33])
bittercress_gdk3b1 <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 2
rda_file <- load(rda_files[34])
bittercress_gdk3b2 <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 3
rda_file <- load(rda_files[35])
bittercress_gdk3b3 <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()

# 
# plot1 <- levelplot(bittercress_gdk0b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "Bittercress gdk 0 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot2 <- levelplot(bittercress_gdk0b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "Bittercress gdk 0 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot3 <- levelplot(bittercress_gdk0b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "Bittercress gdk 0 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot4 <- levelplot(bittercress_gdk1b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "Bittercress gdk 1 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot5 <- levelplot(bittercress_gdk1b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "Bittercress gdk 1 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot6 <- levelplot(bittercress_gdk1b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "Bittercressgdk 1 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot7 <- levelplot(bittercress_gdk3b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "Bittercress gdk 3 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot8 <- levelplot(bittercress_gdk3b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "Bittercress gdk 3 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot9 <- levelplot(bittercress_gdk3b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "Bittercress gdk 3 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# 
# grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol=3, nrow=3)
# 


par(mfrow=c(3,3))

plot(bittercress_gdk0b1)
title(main = c("Bittercress gdk0b1"))

plot(bittercress_gdk0b2)
title(main = c("Bittercress gdk0b2"))

plot(bittercress_gdk0b3)
title(main = c("Bittercress gdk0b3"))

plot(bittercress_gdk1b1)
title(main = c("Bittercress gdk1b1"))

plot(bittercress_gdk1b2)
title(main = c("Bittercress gdk1b2"))

plot(bittercress_gdk1b3)
title(main = c("Bittercress gdk1b3"))

plot(bittercress_gdk3b1)
title(main = c("Bittercress gdk3b1"))

plot(bittercress_gdk3b2)
title(main = c("Bittercress gdk3b2"))

plot(bittercress_gdk3b3)
title(main = c("Bittercress gdk3b3"))



###################################################################




######################################################

setwd('F:/ENMEval Models')

folders <-list.files()

#######################
# bittersweet RESULTS #
#######################

base <- "F:/ENMEval Models/"
fp <- paste(base, folders[2], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#gdk 0 beta 1
rda_file <- load(rda_files[1])
bittersweet_gdk0b1 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 2
rda_file <- load(rda_files[2])
bittersweet_gdk0b2 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 3
rda_file <- load(rda_files[3])
bittersweet_gdk0b3 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 1
rda_file <- load(rda_files[17])
bittersweet_gdk1b1 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 2
rda_file <- load(rda_files[18])
bittersweet_gdk1b2 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 3
rda_file <- load(rda_files[19])
bittersweet_gdk1b3 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 1
rda_file <- load(rda_files[33])
bittersweet_gdk3b1 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 2
rda_file <- load(rda_files[34])
bittersweet_gdk3b2 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 3
rda_file <- load(rda_files[35])
bittersweet_gdk3b3 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()

# 
# plot1 <- levelplot(bittersweet_gdk0b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "bittersweet gdk 0 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot2 <- levelplot(bittersweet_gdk0b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "bittersweet gdk 0 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot3 <- levelplot(bittersweet_gdk0b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "bittersweet gdk 0 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot4 <- levelplot(bittersweet_gdk1b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "bittersweet gdk 1 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot5 <- levelplot(bittersweet_gdk1b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "bittersweet gdk 1 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot6 <- levelplot(bittersweet_gdk1b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "bittersweetgdk 1 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot7 <- levelplot(bittersweet_gdk3b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "bittersweet gdk 3 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot8 <- levelplot(bittersweet_gdk3b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "bittersweet gdk 3 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot9 <- levelplot(bittersweet_gdk3b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "bittersweet gdk 3 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# 
# grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol=3, nrow=3)
# 


par(mfrow=c(3,3))

plot(bittersweet_gdk0b1)
title(main = c("bittersweet gdk0b1"))

plot(bittersweet_gdk0b2)
title(main = c("bittersweet gdk0b2"))

plot(bittersweet_gdk0b3)
title(main = c("bittersweet gdk0b3"))

plot(bittersweet_gdk1b1)
title(main = c("bittersweet gdk1b1"))

plot(bittersweet_gdk1b2)
title(main = c("bittersweet gdk1b2"))

plot(bittersweet_gdk1b3)
title(main = c("bittersweet gdk1b3"))

plot(bittersweet_gdk3b1)
title(main = c("bittersweet gdk3b1"))

plot(bittersweet_gdk3b2)
title(main = c("bittersweet gdk3b2"))

plot(bittersweet_gdk3b3)
title(main = c("bittersweet gdk3b3"))



###################################################################





######################################################

setwd('F:/ENMEval Models')

folders <-list.files()

#######################
# hops RESULTS #
#######################

base <- "F:/ENMEval Models/"
fp <- paste(base, folders[3], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#gdk 0 beta 1
rda_file <- load(rda_files[1])
hops_gdk0b1 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 2
rda_file <- load(rda_files[2])
hops_gdk0b2 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 3
rda_file <- load(rda_files[3])
hops_gdk0b3 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 1
rda_file <- load(rda_files[17])
hops_gdk1b1 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 2
rda_file <- load(rda_files[18])
hops_gdk1b2 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 3
rda_file <- load(rda_files[19])
hops_gdk1b3 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 1
rda_file <- load(rda_files[33])
hops_gdk3b1 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 2
rda_file <- load(rda_files[34])
hops_gdk3b2 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 3
rda_file <- load(rda_files[35])
hops_gdk3b3 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()

# 
# plot1 <- levelplot(hops_gdk0b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "hops gdk 0 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot2 <- levelplot(hops_gdk0b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "hops gdk 0 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot3 <- levelplot(hops_gdk0b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "hops gdk 0 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot4 <- levelplot(hops_gdk1b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "hops gdk 1 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot5 <- levelplot(hops_gdk1b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "hops gdk 1 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot6 <- levelplot(hops_gdk1b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "hopsgdk 1 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot7 <- levelplot(hops_gdk3b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "hops gdk 3 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot8 <- levelplot(hops_gdk3b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "hops gdk 3 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot9 <- levelplot(hops_gdk3b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "hops gdk 3 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# 
# grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol=3, nrow=3)
# 


par(mfrow=c(3,3))

plot(hops_gdk0b1)
title(main = c("hops gdk0b1"))

plot(hops_gdk0b2)
title(main = c("hops gdk0b2"))

plot(hops_gdk0b3)
title(main = c("hops gdk0b3"))

plot(hops_gdk1b1)
title(main = c("hops gdk1b1"))

plot(hops_gdk1b2)
title(main = c("hops gdk1b2"))

plot(hops_gdk1b3)
title(main = c("hops gdk1b3"))

plot(hops_gdk3b1)
title(main = c("hops gdk3b1"))

plot(hops_gdk3b2)
title(main = c("hops gdk3b2"))

plot(hops_gdk3b3)
title(main = c("hops gdk3b3"))



###################################################################






######################################################

setwd('F:/ENMEval Models')

folders <-list.files()

#######################
# knapweed RESULTS #
#######################

base <- "F:/ENMEval Models/"
fp <- paste(base, folders[4], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#gdk 0 beta 1
rda_file <- load(rda_files[1])
knapweed_gdk0b1 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 2
rda_file <- load(rda_files[2])
knapweed_gdk0b2 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 3
rda_file <- load(rda_files[3])
knapweed_gdk0b3 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 1
rda_file <- load(rda_files[17])
knapweed_gdk1b1 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 2
rda_file <- load(rda_files[18])
knapweed_gdk1b2 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 3
rda_file <- load(rda_files[19])
knapweed_gdk1b3 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 1
rda_file <- load(rda_files[33])
knapweed_gdk3b1 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 2
rda_file <- load(rda_files[34])
knapweed_gdk3b2 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 3
rda_file <- load(rda_files[35])
knapweed_gdk3b3 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()

# 
# plot1 <- levelplot(knapweed_gdk0b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "knapweed gdk 0 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot2 <- levelplot(knapweed_gdk0b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "knapweed gdk 0 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot3 <- levelplot(knapweed_gdk0b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "knapweed gdk 0 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot4 <- levelplot(knapweed_gdk1b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "knapweed gdk 1 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot5 <- levelplot(knapweed_gdk1b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "knapweed gdk 1 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot6 <- levelplot(knapweed_gdk1b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "knapweedgdk 1 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot7 <- levelplot(knapweed_gdk3b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "knapweed gdk 3 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot8 <- levelplot(knapweed_gdk3b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "knapweed gdk 3 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot9 <- levelplot(knapweed_gdk3b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "knapweed gdk 3 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# 
# grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol=3, nrow=3)
# 


par(mfrow=c(3,3))

plot(knapweed_gdk0b1)
title(main = c("knapweed gdk0b1"))

plot(knapweed_gdk0b2)
title(main = c("knapweed gdk0b2"))

plot(knapweed_gdk0b3)
title(main = c("knapweed gdk0b3"))

plot(knapweed_gdk1b1)
title(main = c("knapweed gdk1b1"))

plot(knapweed_gdk1b2)
title(main = c("knapweed gdk1b2"))

plot(knapweed_gdk1b3)
title(main = c("knapweed gdk1b3"))

plot(knapweed_gdk3b1)
title(main = c("knapweed gdk3b1"))

plot(knapweed_gdk3b2)
title(main = c("knapweed gdk3b2"))

plot(knapweed_gdk3b3)
title(main = c("knapweed gdk3b3"))



###################################################################






######################################################

setwd('F:/ENMEval Models')

folders <-list.files()

#######################
# parsnip RESULTS #
#######################

base <- "F:/ENMEval Models/"
fp <- paste(base, folders[1], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#gdk 0 beta 1
rda_file <- load(rda_files[1])
parsnip_gdk0b1 <- clusterR(bioclim.parsnip.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 2
rda_file <- load(rda_files[2])
parsnip_gdk0b2 <- clusterR(bioclim.parsnip.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 3
rda_file <- load(rda_files[3])
parsnip_gdk0b3 <- clusterR(bioclim.parsnip.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 1
rda_file <- load(rda_files[17])
parsnip_gdk1b1 <- clusterR(bioclim.parsnip.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 2
rda_file <- load(rda_files[18])
parsnip_gdk1b2 <- clusterR(bioclim.parsnip.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 3
rda_file <- load(rda_files[19])
parsnip_gdk1b3 <- clusterR(bioclim.parsnip.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 1
rda_file <- load(rda_files[33])
parsnip_gdk3b1 <- clusterR(bioclim.parsnip.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 2
rda_file <- load(rda_files[34])
parsnip_gdk3b2 <- clusterR(bioclim.parsnip.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 3
rda_file <- load(rda_files[35])
parsnip_gdk3b3 <- clusterR(bioclim.parsnip.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()

# 
# plot1 <- levelplot(parsnip_gdk0b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "parsnip gdk 0 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot2 <- levelplot(parsnip_gdk0b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "parsnip gdk 0 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot3 <- levelplot(parsnip_gdk0b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "parsnip gdk 0 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot4 <- levelplot(parsnip_gdk1b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "parsnip gdk 1 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot5 <- levelplot(parsnip_gdk1b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "parsnip gdk 1 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot6 <- levelplot(parsnip_gdk1b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "parsnipgdk 1 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot7 <- levelplot(parsnip_gdk3b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "parsnip gdk 3 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot8 <- levelplot(parsnip_gdk3b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "parsnip gdk 3 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot9 <- levelplot(parsnip_gdk3b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "parsnip gdk 3 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# 
# grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol=3, nrow=3)
# 


par(mfrow=c(3,3))

plot(parsnip_gdk0b1)
title(main = c("parsnip gdk0b1"))

plot(parsnip_gdk0b2)
title(main = c("parsnip gdk0b2"))

plot(parsnip_gdk0b3)
title(main = c("parsnip gdk0b3"))

plot(parsnip_gdk1b1)
title(main = c("parsnip gdk1b1"))

plot(parsnip_gdk1b2)
title(main = c("parsnip gdk1b2"))

plot(parsnip_gdk1b3)
title(main = c("parsnip gdk1b3"))

plot(parsnip_gdk3b1)
title(main = c("parsnip gdk3b1"))

plot(parsnip_gdk3b2)
title(main = c("parsnip gdk3b2"))

plot(parsnip_gdk3b3)
title(main = c("parsnip gdk3b3"))



###################################################################




######################################################

setwd('F:/ENMEval Models')

folders <-list.files()

#######################
# swallowwort RESULTS #
#######################

base <- "F:/ENMEval Models/"
fp <- paste(base, folders[2], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#gdk 0 beta 1
rda_file <- load(rda_files[1])
swallowwort_gdk0b1 <- clusterR(bioclim.sw.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 2
rda_file <- load(rda_files[2])
swallowwort_gdk0b2 <- clusterR(bioclim.sw.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 3
rda_file <- load(rda_files[3])
swallowwort_gdk0b3 <- clusterR(bioclim.sw.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 1
rda_file <- load(rda_files[17])
swallowwort_gdk1b1 <- clusterR(bioclim.sw.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 2
rda_file <- load(rda_files[18])
swallowwort_gdk1b2 <- clusterR(bioclim.sw.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 3
rda_file <- load(rda_files[19])
swallowwort_gdk1b3 <- clusterR(bioclim.sw.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 1
rda_file <- load(rda_files[33])
swallowwort_gdk3b1 <- clusterR(bioclim.sw.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 2
rda_file <- load(rda_files[34])
swallowwort_gdk3b2 <- clusterR(bioclim.sw.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 3
rda_file <- load(rda_files[35])
swallowwort_gdk3b3 <- clusterR(bioclim.sw.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()

# 
# plot1 <- levelplot(swallowwort_gdk0b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "swallowwort gdk 0 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot2 <- levelplot(swallowwort_gdk0b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "swallowwort gdk 0 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot3 <- levelplot(swallowwort_gdk0b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "swallowwort gdk 0 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot4 <- levelplot(swallowwort_gdk1b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "swallowwort gdk 1 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot5 <- levelplot(swallowwort_gdk1b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "swallowwort gdk 1 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot6 <- levelplot(swallowwort_gdk1b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "swallowwortgdk 1 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot7 <- levelplot(swallowwort_gdk3b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "swallowwort gdk 3 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot8 <- levelplot(swallowwort_gdk3b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "swallowwort gdk 3 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot9 <- levelplot(swallowwort_gdk3b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "swallowwort gdk 3 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# 
# grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol=3, nrow=3)
# 


par(mfrow=c(3,3))

plot(swallowwort_gdk0b1)
title(main = c("swallowwort gdk0b1"))

plot(swallowwort_gdk0b2)
title(main = c("swallowwort gdk0b2"))

plot(swallowwort_gdk0b3)
title(main = c("swallowwort gdk0b3"))

plot(swallowwort_gdk1b1)
title(main = c("swallowwort gdk1b1"))

plot(swallowwort_gdk1b2)
title(main = c("swallowwort gdk1b2"))

plot(swallowwort_gdk1b3)
title(main = c("swallowwort gdk1b3"))

plot(swallowwort_gdk3b1)
title(main = c("swallowwort gdk3b1"))

plot(swallowwort_gdk3b2)
title(main = c("swallowwort gdk3b2"))

plot(swallowwort_gdk3b3)
title(main = c("swallowwort gdk3b3"))



###################################################################






######################################################

setwd('F:/ENMEval Models')

folders <-list.files()

#######################
# tansy RESULTS #
#######################

base <- "F:/ENMEval Models/"
fp <- paste(base, folders[3], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#gdk 0 beta 1
rda_file <- load(rda_files[1])
tansy_gdk0b1 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 2
rda_file <- load(rda_files[2])
tansy_gdk0b2 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 3
rda_file <- load(rda_files[3])
tansy_gdk0b3 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 1
rda_file <- load(rda_files[17])
tansy_gdk1b1 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 2
rda_file <- load(rda_files[18])
tansy_gdk1b2 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 3
rda_file <- load(rda_files[19])
tansy_gdk1b3 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 1
rda_file <- load(rda_files[33])
tansy_gdk3b1 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 2
rda_file <- load(rda_files[34])
tansy_gdk3b2 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 3
rda_file <- load(rda_files[35])
tansy_gdk3b3 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()

# 
# plot1 <- levelplot(tansy_gdk0b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "tansy gdk 0 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot2 <- levelplot(tansy_gdk0b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "tansy gdk 0 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot3 <- levelplot(tansy_gdk0b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "tansy gdk 0 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot4 <- levelplot(tansy_gdk1b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "tansy gdk 1 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot5 <- levelplot(tansy_gdk1b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "tansy gdk 1 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot6 <- levelplot(tansy_gdk1b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "tansygdk 1 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot7 <- levelplot(tansy_gdk3b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "tansy gdk 3 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot8 <- levelplot(tansy_gdk3b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "tansy gdk 3 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot9 <- levelplot(tansy_gdk3b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "tansy gdk 3 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# 
# grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol=3, nrow=3)
# 


par(mfrow=c(3,3))

plot(tansy_gdk0b1)
title(main = c("tansy gdk0b1"))

plot(tansy_gdk0b2)
title(main = c("tansy gdk0b2"))

plot(tansy_gdk0b3)
title(main = c("tansy gdk0b3"))

plot(tansy_gdk1b1)
title(main = c("tansy gdk1b1"))

plot(tansy_gdk1b2)
title(main = c("tansy gdk1b2"))

plot(tansy_gdk1b3)
title(main = c("tansy gdk1b3"))

plot(tansy_gdk3b1)
title(main = c("tansy gdk3b1"))

plot(tansy_gdk3b2)
title(main = c("tansy gdk3b2"))

plot(tansy_gdk3b3)
title(main = c("tansy gdk3b3"))



###################################################################






######################################################

setwd('F:/ENMEval Models')

folders <-list.files()

#######################
# teasel RESULTS #
#######################

base <- "F:/ENMEval Models/"
fp <- paste(base, folders[8], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#gdk 0 beta 1
rda_file <- load(rda_files[1])
teasel_gdk0b1 <- clusterR(bioclim.ct.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 2
rda_file <- load(rda_files[2])
teasel_gdk0b2 <- clusterR(bioclim.ct.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 3
rda_file <- load(rda_files[3])
teasel_gdk0b3 <- clusterR(bioclim.ct.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 1
rda_file <- load(rda_files[17])
teasel_gdk1b1 <- clusterR(bioclim.ct.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 2
rda_file <- load(rda_files[18])
teasel_gdk1b2 <- clusterR(bioclim.ct.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 3
rda_file <- load(rda_files[19])
teasel_gdk1b3 <- clusterR(bioclim.ct.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 1
rda_file <- load(rda_files[33])
teasel_gdk3b1 <- clusterR(bioclim.ct.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 2
rda_file <- load(rda_files[34])
teasel_gdk3b2 <- clusterR(bioclim.ct.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 3
rda_file <- load(rda_files[35])
teasel_gdk3b3 <- clusterR(bioclim.ct.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()

# 
# plot1 <- levelplot(teasel_gdk0b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "teasel gdk 0 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot2 <- levelplot(teasel_gdk0b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "teasel gdk 0 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot3 <- levelplot(teasel_gdk0b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "teasel gdk 0 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot4 <- levelplot(teasel_gdk1b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "teasel gdk 1 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot5 <- levelplot(teasel_gdk1b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "teasel gdk 1 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot6 <- levelplot(teasel_gdk1b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "teaselgdk 1 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot7 <- levelplot(teasel_gdk3b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "teasel gdk 3 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot8 <- levelplot(teasel_gdk3b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "teasel gdk 3 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot9 <- levelplot(teasel_gdk3b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "teasel gdk 3 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# 
# grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol=3, nrow=3)
# 


par(mfrow=c(3,3))

plot(teasel_gdk0b1)
title(main = c("teasel gdk0b1"))

plot(teasel_gdk0b2)
title(main = c("teasel gdk0b2"))

plot(teasel_gdk0b3)
title(main = c("teasel gdk0b3"))

plot(teasel_gdk1b1)
title(main = c("teasel gdk1b1"))

plot(teasel_gdk1b2)
title(main = c("teasel gdk1b2"))

plot(teasel_gdk1b3)
title(main = c("teasel gdk1b3"))

plot(teasel_gdk3b1)
title(main = c("teasel gdk3b1"))

plot(teasel_gdk3b2)
title(main = c("teasel gdk3b2"))

plot(teasel_gdk3b3)
title(main = c("teasel gdk3b3"))



###################################################################






######################################################

setwd('F:/ENMEval Models')

folders <-list.files()

#######################
# toadflax RESULTS #
#######################

base <- "F:/ENMEval Models/"
fp <- paste(base, folders[9], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#gdk 0 beta 1
rda_file <- load(rda_files[1])
toadflax_gdk0b1 <- clusterR(bioclim.tf.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 2
rda_file <- load(rda_files[2])
toadflax_gdk0b2 <- clusterR(bioclim.tf.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 beta 3
rda_file <- load(rda_files[3])
toadflax_gdk0b3 <- clusterR(bioclim.tf.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 1
rda_file <- load(rda_files[17])
toadflax_gdk1b1 <- clusterR(bioclim.tf.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 2
rda_file <- load(rda_files[18])
toadflax_gdk1b2 <- clusterR(bioclim.tf.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 beta 3
rda_file <- load(rda_files[19])
toadflax_gdk1b3 <- clusterR(bioclim.tf.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 1
rda_file <- load(rda_files[33])
toadflax_gdk3b1 <- clusterR(bioclim.tf.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 2
rda_file <- load(rda_files[34])
toadflax_gdk3b2 <- clusterR(bioclim.tf.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 3 beta 3
rda_file <- load(rda_files[35])
toadflax_gdk3b3 <- clusterR(bioclim.tf.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()

# 
# plot1 <- levelplot(toadflax_gdk0b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "toadflax gdk 0 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot2 <- levelplot(toadflax_gdk0b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "toadflax gdk 0 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot3 <- levelplot(toadflax_gdk0b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "toadflax gdk 0 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot4 <- levelplot(toadflax_gdk1b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "toadflax gdk 1 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot5 <- levelplot(toadflax_gdk1b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "toadflax gdk 1 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot6 <- levelplot(toadflax_gdk1b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "toadflaxgdk 1 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot7 <- levelplot(toadflax_gdk3b1, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "toadflax gdk 3 b 1") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot8 <- levelplot(toadflax_gdk3b2, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "toadflax gdk 3 b 2") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# plot9 <- levelplot(toadflax_gdk3b3, margin = FALSE, maxpixels = 1e5, col.regions = magma(256, begin = 0, end = 1), main = "toadflax gdk 3 b 3") + spplot(states.lines.crop, col.regions = "slateblue1") #+ spplot(occurrence_current, pch = 19, cex = 0.05, col.regions = "red")
# 
# 
# grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol=3, nrow=3)
# 


par(mfrow=c(3,3))

plot(toadflax_gdk0b1)
title(main = c("toadflax gdk0b1"))

plot(toadflax_gdk0b2)
title(main = c("toadflax gdk0b2"))

plot(toadflax_gdk0b3)
title(main = c("toadflax gdk0b3"))

plot(toadflax_gdk1b1)
title(main = c("toadflax gdk1b1"))

plot(toadflax_gdk1b2)
title(main = c("toadflax gdk1b2"))

plot(toadflax_gdk1b3)
title(main = c("toadflax gdk1b3"))

plot(toadflax_gdk3b1)
title(main = c("toadflax gdk3b1"))

plot(toadflax_gdk3b2)
title(main = c("toadflax gdk3b2"))

plot(toadflax_gdk3b3)
title(main = c("toadflax gdk3b3"))



###################################################################











































#######################################################################


##### RESULTS May 20 2019 12 Panel Figures ###############


#######################################################################


setwd('F:/ENMEval Models')

folders <-list.files()

#######################
# bittercress RESULTS #
#######################

base <- "F:/ENMEval Models/"
fp <- paste(base, folders[1], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#block models

#### Model Group 1 - overall max model####

#max auc model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_0_ds_1km_method_block_model_2.rda")
max_auc <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max tss model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_0_ds_1km_method_block_model_2.rda")
max_tss <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max boyce model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_0_ds_10km_method_block_model_3.rda")
max_boyce <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max sensitivity model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_1_ds_10km_method_block_model_4.rda")
max_sens <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)


#### Model Group 2 - max boi models (block only)####

#max auc BOI
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_0_ds_1km_method_block_model_1.rda")
max_auc_boi <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max tss BOI
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_0_ds_1km_method_block_model_4.rda")
max_tss_boi <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max boyce BOI
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_0_ds_10km_method_block_model_3.rda")
max_boyce_boi <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max sensitivity BOI
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_1_ds_10km_method_block_model_4.rda")
max_sens_boi <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)


#### Model Group 3 - min variation models ####

#max auc model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_3_ds_10km_method_block_model_1.rda")
min_auc_var <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max tss model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_3_ds_1km_method_block_model_1.rda")
min_tss_var <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max boyce model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_0_ds_10km_method_block_model_1.rda")
min_boyce_var <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max sensitivity model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_3_ds_1km_method_block_model_2.rda")
min_sens_var <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()


#bittercress

bias_list <- get(paste(species.names[8], "_bias_list", sep="")) #get the species_bias_list object 

#prepare species occurrence and background data
d <- prepareData(species.file.paths[8], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

blocks <- get.block(d[[1]], d[[2]])

occ_blocks <- cbind(d[[1]], blocks$occ.grp)

block1 <- occ_blocks[which(occ_blocks[,3]==1),]

block2 <- occ_blocks[which(occ_blocks[,3]==2),]

block3 <- occ_blocks[which(occ_blocks[,3]==3),]

block4 <- occ_blocks[which(occ_blocks[,3]==4),]

#start plotting

par(mfrow=c(3,5))

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Bittercress Block")

points(block4, col="grey", pch=20)
points(block3, col="purple", pch=20)
points(block2, col="green", pch=20)
points(block1, col="red", pch=20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))

plot(max_auc)
title(main=c("Bittercress Max Overall Model AUC"))
legend("topright", legend=c("GDK 0 Downsample 1 Beta 2"))
points(bi)

plot(max_tss)
title(main=c("Bittercress Max Overall Model TSS"))
legend("topright", legend=c("GDK 0 Downsample 1 Beta 2"))

plot(max_boyce)
title(main=c("Bittercress Max Overall Model Boyce"))
legend("topright", legend=c("GDK 0 Downsample 10 Beta 3"))

plot(max_sens)
title(main=c("Bittercress Max Overall Model Sensitivity"))
legend("topright", legend=c("GDK 1 Downsample 10 Beta 4"))


plot(mean.temp.bio1, col="gray", legend=FALSE, main="Bittercress Presences")
points(d[[1]], col = "black", pch =20)

#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))


plot(max_auc_boi)
title(main=c("Bittercress Max AUC Block 3"))
legend("topright", legend=c("GDK 0 Downsample 1 Beta 1"))

plot(max_tss_boi)
title(main=c("Bittercress Max TSS Block 3"))
legend("topright", legend=c("GDK 0 Downsample 1 Beta 4"))

plot(max_boyce_boi)
title(main=c("Bittercress Max Boyce Block 3"))
legend("topright", legend=c("GDK 0 Downsample 10 Beta 3"))

plot(max_sens_boi)
title(main=c("Bittercress Max Sensitivity Block 3"))
legend("topright", legend=c("GDK 1 Downsample 10 Beta 4"))

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Bittercress Background")
points(d[[2]], col = "black", pch = 20)


#points(d[[1]], pch=21, bg=blocks$occ.grp)

legend("topright", legend=c("block1", "block2", "block3", "block4"), pch=20, col=c("red", "green", "purple", "black"))


plot(min_auc_var)
title(main=c("Bittercress Min Variation AUC Blocks"))
legend("topright", legend=c("GDK 3 Downsample 10 Beta 1"))

plot(min_tss_var)
title(main=c("Bittercress Min Variaton TSS Blocks"))
legend("topright", legend=c("GDK 3 Downsample 1 Beta 1"))

plot(min_boyce_var)
title(main=c("Bittercress Min Variation Boyce Blocks"))
legend("topright", legend=c("GDK 0 Downsample 10 Beta 1"))

plot(min_sens_var)
title(main=c("Bittercress Min Variation Sensitivity Blocks"))
legend("topright", legend=c("GDK 3 Downsample 1 Beta 2"))




setwd('F:/ENMEval Models')

folders <-list.files()

#######################
# bittercress RESULTS #
#######################

base <- "F:/ENMEval Models/"
fp <- paste(base, folders[1], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#checkerboard models

#### Model Group 1 - overall max model####

#max auc model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_0_ds_1km_method_checkerboard2_model_1.rda")
max_auc <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max tss model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_0_ds_1km_method_checkerboard2_model_1.rda")
max_tss <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max boyce model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_1_ds_1km_method_checkerboard2_model_2.rda")
max_boyce <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max sensitivity model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_0_ds_1km_method_checkerboard2_model_1.rda")
max_sens <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)


#### Model Group 2 - min variation models ####

#max auc model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_0_ds_10km_method_checkerboard2_model_2.rda")
min_auc_var <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max tss model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_0_ds_10km_method_checkerboard2_model_1.rda")
min_tss_var <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max boyce model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_3_ds_10km_method_checkerboard2_model_3.rda")
min_boyce_var <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

#max sensitivity model
rda_file <- load("F:/ENMEval Models/bittercress/species_bittercress_gdk_3_ds_10km_method_checkerboard2_model_4.rda")
min_sens_var <- clusterR(bioclim.nb.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()


par(mfrow=c(2, 5))


plot(mean.temp.bio1, col="gray", legend=FALSE, main="Bittercress Presences")
points(d[[1]], col = "black", pch =20)

plot(max_auc)
title(main=c("Bittercress Max Overall Model AUC"))
legend("topright", legend=c("GDK 0 Downsample 1 Beta 1"))
points(bi)

plot(max_tss)
title(main=c("Bittercress Max Overall Model TSS"))
legend("topright", legend=c("GDK 0 Downsample 1 Beta 1"))

plot(max_boyce)
title(main=c("Bittercress Max Overall Model Boyce"))
legend("topright", legend=c("GDK 1 Downsample 1 Beta 2"))

plot(max_sens)
title(main=c("Bittercress Max Overall Model Sensitivity"))
legend("topright", legend=c("GDK 1 Downsample 1 Beta 1"))

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Bittercress Background")
points(d[[2]], col = "black", pch = 20)

plot(min_auc_var)
title(main=c("Bittercress Min Variation AUC Blocks"))
legend("topright", legend=c("GDK 0 Downsample 10 Beta 2"))

plot(min_tss_var)
title(main=c("Bittercress Min Variaton TSS Blocks"))
legend("topright", legend=c("GDK 0 Downsample 10 Beta 1"))

plot(min_boyce_var)
title(main=c("Bittercress Min Variation Boyce Blocks"))
legend("topright", legend=c("GDK 3 Downsample 10 Beta 3"))

plot(min_sens_var)
title(main=c("Bittercress Min Variation Sensitivity Blocks"))
legend("topright", legend=c("GDK 3 Downsample 10 Beta 4"))




###########################################################################################################


















#################################################################################################################


#################################
 
##### Threshold Area Plots #####

#################################

##### June 12 2019 #####

#Testing with knapweed


setwd('F:/ENMEval Models')

folders <-list.files()


base <- "F:/ENMEval Models/"
fp <- paste(base, folders[4], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#block

#gdk 0 ds 10
rda_file <- load(rda_files[1])
knapweed_gdk0_ds10 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 ds 10
rda_file <- load(rda_files[17])
knapweed_gdk1_ds10 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)


#gdk 3 ds 10
rda_file <- load(rda_files[33])
knapweed_gdk3_ds10 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 ds 1
rda_file <- load(rda_files[9])
knapweed_gdk0_ds1 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 ds 1
rda_file <- load(rda_files[25])
knapweed_gdk1_ds1 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)


#gdk 3 ds 1
rda_file <- load(rda_files[41])
knapweed_gdk3_ds1 <- clusterR(bioclim.bk.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()

#test threshold area plot

#gdk0ds10

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- knapweed_gdk0_ds10 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk0ds10 <- as.data.frame(all_threshold_areas)
gdk0ds10 <- rbind(gdk0ds10, rep(1:100))
gdk0ds10 <- t(gdk0ds10)
colnames(gdk0ds10) <- c("area", "threshold_percent")

#gdk1ds10

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- knapweed_gdk1_ds10 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk1ds10 <- as.data.frame(all_threshold_areas)
gdk1ds10 <- rbind(gdk1ds10, rep(1:100))
gdk1ds10 <- t(gdk1ds10)
colnames(gdk1ds10) <- c("area", "threshold_percent")

#gdk3ds10

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- knapweed_gdk3_ds10 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk3ds10 <- as.data.frame(all_threshold_areas)
gdk3ds10 <- rbind(gdk3ds10, rep(1:100))
gdk3ds10 <- t(gdk3ds10)
colnames(gdk3ds10) <- c("area", "threshold_percent")

#gdk0ds1

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- knapweed_gdk0_ds1 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk0ds1 <- as.data.frame(all_threshold_areas)
gdk0ds1 <- rbind(gdk0ds1, rep(1:100))
gdk0ds1 <- t(gdk0ds1)
colnames(gdk0ds1) <- c("area", "threshold_percent")

#gdk1ds1

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- knapweed_gdk1_ds1 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk1ds1 <- as.data.frame(all_threshold_areas)
gdk1ds1 <- rbind(gdk1ds1, rep(1:100))
gdk1ds1 <- t(gdk1ds1)
colnames(gdk1ds1) <- c("area", "threshold_percent")

#gdk3ds1

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- knapweed_gdk3_ds1 = j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk3ds1 <- as.data.frame(all_threshold_areas)
gdk3ds1 <- rbind(gdk3ds1, rep(1:100))
gdk3ds1 <- t(gdk3ds1)
colnames(gdk3ds1) <- c("area", "threshold_percent")


#plotting

par(mfrow=c(4,2))
plot(knapweed_gdk0_ds10, main = c("GDK0 DS 10"))
plot(knapweed_gdk0_ds1, main = c("GDK0 DS 1"))

plot(knapweed_gdk1_ds10, main = c("GDK1 DS 10"))
plot(knapweed_gdk1_ds1, main = c("GDK1 DS 1"))

plot(knapweed_gdk3_ds10, main = c("GDK3 DS 10"))
plot(knapweed_gdk3_ds1, main = c("GDK3 DS 1"))

par(mar=c(5,4,4,4))

#gdk3
plot(gdk3ds10[,2] ~ gdk3ds10[,1], ylab = "Threshold", xlab = "Total pixels", main = c("Threshold Area Plot"), type="l", col = "black", lty = 1)
lines(gdk3ds1[,2] ~ gdk3ds1[,1], type="l", col = "black", lty = 5)

#gdk1
lines(gdk1ds10[,2]~ gdk1ds10[,1], type = "l", col='blue', lty = 1)
lines(gdk1ds1[,2] ~ gdk1ds1[,1], type="l", col = "blue", lty = 5)

#gdk0
lines(gdk0ds10[,2]~ gdk0ds10[,1], type = "l", col='red', lty = 1)
lines(gdk0ds1[,2] ~ gdk0ds1[,1], type="l", col = "red", lty = 5)


legend("topright", legend=c("GDK3 DS10", "GDK3 DS1", "GDK1 D 10", "GDK1 DS1", "GDK0 DS10", "GDK0 DS1"), col = c("black", "black", "blue", "blue", "red", "red"), lty = c(1,5,1,5,1,5))

#knapweed

bias_list <- get(paste(species.names[2], "_bias_list", sep="")) #get the species_bias_list object 

#prepare species occurrence and background data
d <- prepareData(species.file.paths[2], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

plot(mean.temp.bio1, col="gray", legend=FALSE, main="knapweed Occurrences")
points(d[[1]], pch = 19)






#################################################################################################################


#################################

##### Threshold Area Plots #####

#################################

##### June 12 2019 #####

#Testing with bittersweet


setwd('F:/ENMEval Models')

folders <-list.files()


base <- "F:/ENMEval Models/"
fp <- paste(base, folders[2], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#block

#gdk 0 ds 10
rda_file <- load(rda_files[1])
bittersweet_gdk0_ds10 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 ds 10
rda_file <- load(rda_files[17])
bittersweet_gdk1_ds10 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)


#gdk 3 ds 10
rda_file <- load(rda_files[33])
bittersweet_gdk3_ds10 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 ds 1
rda_file <- load(rda_files[9])
bittersweet_gdk0_ds1 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 ds 1
rda_file <- load(rda_files[25])
bittersweet_gdk1_ds1 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)


#gdk 3 ds 1
rda_file <- load(rda_files[41])
bittersweet_gdk3_ds1 <- clusterR(bioclim.ob.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()

#test threshold area plot

#gdk0ds10

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- bittersweet_gdk0_ds10 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk0ds10 <- as.data.frame(all_threshold_areas)
gdk0ds10 <- rbind(gdk0ds10, rep(1:100))
gdk0ds10 <- t(gdk0ds10)
colnames(gdk0ds10) <- c("area", "threshold_percent")

#gdk1ds10

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- bittersweet_gdk1_ds10 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk1ds10 <- as.data.frame(all_threshold_areas)
gdk1ds10 <- rbind(gdk1ds10, rep(1:100))
gdk1ds10 <- t(gdk1ds10)
colnames(gdk1ds10) <- c("area", "threshold_percent")

#gdk3ds10

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- bittersweet_gdk3_ds10 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk3ds10 <- as.data.frame(all_threshold_areas)
gdk3ds10 <- rbind(gdk3ds10, rep(1:100))
gdk3ds10 <- t(gdk3ds10)
colnames(gdk3ds10) <- c("area", "threshold_percent")

#gdk0ds1

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- bittersweet_gdk0_ds1 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk0ds1 <- as.data.frame(all_threshold_areas)
gdk0ds1 <- rbind(gdk0ds1, rep(1:100))
gdk0ds1 <- t(gdk0ds1)
colnames(gdk0ds1) <- c("area", "threshold_percent")

#gdk1ds1

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- bittersweet_gdk1_ds1 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk1ds1 <- as.data.frame(all_threshold_areas)
gdk1ds1 <- rbind(gdk1ds1, rep(1:100))
gdk1ds1 <- t(gdk1ds1)
colnames(gdk1ds1) <- c("area", "threshold_percent")

#gdk3ds1

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- bittersweet_gdk3_ds1 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk3ds1 <- as.data.frame(all_threshold_areas)
gdk3ds1 <- rbind(gdk3ds1, rep(1:100))
gdk3ds1 <- t(gdk3ds1)
colnames(gdk3ds1) <- c("area", "threshold_percent")


#plotting

par(mar=c(5,4,4,4))

par(mfrow=c(4,2))
plot(bittersweet_gdk0_ds10, main = c("GDK0 DS 10"))
plot(bittersweet_gdk0_ds1, main = c("GDK0 DS 1"))

plot(bittersweet_gdk1_ds10, main = c("GDK1 DS 10"))
plot(bittersweet_gdk1_ds1, main = c("GDK1 DS 1"))

plot(bittersweet_gdk3_ds10, main = c("GDK3 DS 10"))
plot(bittersweet_gdk3_ds1, main = c("GDK3 DS 1"))



#gdk3
plot(gdk3ds10[,2] ~ gdk3ds10[,1], ylab = "Threshold", xlab = "Total pixels", main = c("Threshold Area Plot"), type="l", col = "black", lty = 1, xlim = c(0, 3e+07))
lines(gdk3ds1[,2] ~ gdk3ds1[,1], type="l", col = "black", lty = 5)

#gdk1
lines(gdk1ds10[,2]~ gdk1ds10[,1], type = "l", col='blue', lty = 1)
lines(gdk1ds1[,2] ~ gdk1ds1[,1], type="l", col = "blue", lty = 5)

#gdk0
lines(gdk0ds10[,2]~ gdk0ds10[,1], type = "l", col='red', lty = 1)
lines(gdk0ds1[,2] ~ gdk0ds1[,1], type="l", col = "red", lty = 5)


legend("topright", legend=c("GDK3 DS10", "GDK3 DS1", "GDK1 DS10", "GDK1 DS1", "GDK0 DS10", "GDK0 DS1"), col = c("black", "black", "blue", "blue", "red", "red"), lty = c(1,5,1,5,1,5))

#bittersweet

bias_list <- get(paste(species.names[9], "_bias_list", sep="")) #get the species_bias_list object 

#prepare species occurrence and background data
d <- prepareData(species.file.paths[9], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

plot(mean.temp.bio1, col="gray", legend=FALSE, main="Bittersweet Occurrences")
points(d[[1]], pch = 19)




#################################################################################################################


#################################

##### Threshold Area Plots #####

#################################

##### June 12 2019 #####

#Testing with bittersweet


setwd('F:/ENMEval Models')

folders <-list.files()


base <- "F:/ENMEval Models/"
fp <- paste(base, folders[9], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#block

#gdk 0 ds 10
rda_file <- load(rda_files[1])
tansy_gdk0_ds10 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 ds 10
rda_file <- load(rda_files[17])
tansy_gdk1_ds10 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)


#gdk 3 ds 10
rda_file <- load(rda_files[33])
tansy_gdk3_ds10 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 ds 1
rda_file <- load(rda_files[9])
tansy_gdk0_ds1 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 ds 1
rda_file <- load(rda_files[25])
tansy_gdk1_ds1 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)


#gdk 3 ds 1
rda_file <- load(rda_files[41])
tansy_gdk3_ds1 <- clusterR(bioclim.tansy.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()

#test threshold area plot

#gdk0ds10

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- tansy_gdk0_ds10 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk0ds10 <- as.data.frame(all_threshold_areas)
gdk0ds10 <- rbind(gdk0ds10, rep(1:100))
gdk0ds10 <- t(gdk0ds10)
colnames(gdk0ds10) <- c("area", "threshold_percent")

#gdk1ds10

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- tansy_gdk1_ds10 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk1ds10 <- as.data.frame(all_threshold_areas)
gdk1ds10 <- rbind(gdk1ds10, rep(1:100))
gdk1ds10 <- t(gdk1ds10)
colnames(gdk1ds10) <- c("area", "threshold_percent")

#gdk3ds10

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- tansy_gdk3_ds10 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk3ds10 <- as.data.frame(all_threshold_areas)
gdk3ds10 <- rbind(gdk3ds10, rep(1:100))
gdk3ds10 <- t(gdk3ds10)
colnames(gdk3ds10) <- c("area", "threshold_percent")

#gdk0ds1

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- tansy_gdk0_ds1 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk0ds1 <- as.data.frame(all_threshold_areas)
gdk0ds1 <- rbind(gdk0ds1, rep(1:100))
gdk0ds1 <- t(gdk0ds1)
colnames(gdk0ds1) <- c("area", "threshold_percent")

#gdk1ds1

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- tansy_gdk1_ds1 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk1ds1 <- as.data.frame(all_threshold_areas)
gdk1ds1 <- rbind(gdk1ds1, rep(1:100))
gdk1ds1 <- t(gdk1ds1)
colnames(gdk1ds1) <- c("area", "threshold_percent")

#gdk3ds1

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- tansy_gdk3_ds1 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk3ds1 <- as.data.frame(all_threshold_areas)
gdk3ds1 <- rbind(gdk3ds1, rep(1:100))
gdk3ds1 <- t(gdk3ds1)
colnames(gdk3ds1) <- c("area", "threshold_percent")


#plotting

par(mar=c(5,4,4,4))

par(mfrow=c(4,2))
plot(tansy_gdk0_ds10, main = c("GDK0 DS 10"))
plot(tansy_gdk0_ds1, main = c("GDK0 DS 1"))

plot(tansy_gdk1_ds10, main = c("GDK1 DS 10"))
plot(tansy_gdk1_ds1, main = c("GDK1 DS 1"))

plot(tansy_gdk3_ds10, main = c("GDK3 DS 10"))
plot(tansy_gdk3_ds1, main = c("GDK3 DS 1"))



#gdk3
plot(gdk3ds10[,2] ~ gdk3ds10[,1], ylab = "Threshold", xlab = "Total pixels", main = c("Threshold Area Plot"), type="l", col = "black", lty = 1, xlim = c(0, 3e+07))
lines(gdk3ds1[,2] ~ gdk3ds1[,1], type="l", col = "black", lty = 5)

#gdk1
lines(gdk1ds10[,2]~ gdk1ds10[,1], type = "l", col='blue', lty = 1)
lines(gdk1ds1[,2] ~ gdk1ds1[,1], type="l", col = "blue", lty = 5)

#gdk0
lines(gdk0ds10[,2]~ gdk0ds10[,1], type = "l", col='red', lty = 1)
lines(gdk0ds1[,2] ~ gdk0ds1[,1], type="l", col = "red", lty = 5)


legend("topright", legend=c("GDK3 DS10", "GDK3 DS1", "GDK1 DS10", "GDK1 DS1", "GDK0 DS10", "GDK0 DS1"), col = c("black", "black", "blue", "blue", "red", "red"), lty = c(1,5,1,5,1,5))

#tansy

bias_list <- get(paste(species.names[6], "_bias_list", sep="")) #get the species_bias_list object 

#prepare species occurrence and background data
d <- prepareData(species.file.paths[6], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

plot(mean.temp.bio1, col="gray", legend=FALSE, main="tansy Occurrences")
points(d[[1]], pch = 19)








#################################################################################################################


#################################

##### Threshold Area Plots #####

#################################

##### June 12 2019 #####

#Testing with bittersweet


setwd('F:/ENMEval Models')

folders <-list.files()


base <- "F:/ENMEval Models/"
fp <- paste(base, folders[4], sep="")
rda_files <- list.files(fp, full.names = TRUE)


beginCluster()

#block

#gdk 0 ds 10
rda_file <- load(rda_files[1])
hops_gdk0_ds10 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 ds 10
rda_file <- load(rda_files[17])
hops_gdk1_ds10 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)


#gdk 3 ds 10
rda_file <- load(rda_files[33])
hops_gdk3_ds10 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 0 ds 1
rda_file <- load(rda_files[9])
hops_gdk0_ds1 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)

#gdk 1 ds 1
rda_file <- load(rda_files[25])
hops_gdk1_ds1 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)


#gdk 3 ds 1
rda_file <- load(rda_files[41])
hops_gdk3_ds1 <- clusterR(bioclim.jh.3, raster::predict, args = list(model = best.model))
rm(best.model)

endCluster()

#test threshold area plot

#gdk0ds10

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- hops_gdk0_ds10 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk0ds10 <- as.data.frame(all_threshold_areas)
gdk0ds10 <- rbind(gdk0ds10, rep(1:100))
gdk0ds10 <- t(gdk0ds10)
colnames(gdk0ds10) <- c("area", "threshold_percent")

#gdk1ds10

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- hops_gdk1_ds10 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk1ds10 <- as.data.frame(all_threshold_areas)
gdk1ds10 <- rbind(gdk1ds10, rep(1:100))
gdk1ds10 <- t(gdk1ds10)
colnames(gdk1ds10) <- c("area", "threshold_percent")

#gdk3ds10

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- hops_gdk3_ds10 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk3ds10 <- as.data.frame(all_threshold_areas)
gdk3ds10 <- rbind(gdk3ds10, rep(1:100))
gdk3ds10 <- t(gdk3ds10)
colnames(gdk3ds10) <- c("area", "threshold_percent")

#gdk0ds1

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- hops_gdk0_ds1 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk0ds1 <- as.data.frame(all_threshold_areas)
gdk0ds1 <- rbind(gdk0ds1, rep(1:100))
gdk0ds1 <- t(gdk0ds1)
colnames(gdk0ds1) <- c("area", "threshold_percent")

#gdk1ds1

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- hops_gdk1_ds1 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk1ds1 <- as.data.frame(all_threshold_areas)
gdk1ds1 <- rbind(gdk1ds1, rep(1:100))
gdk1ds1 <- t(gdk1ds1)
colnames(gdk1ds1) <- c("area", "threshold_percent")

#gdk3ds1

all_threshold_areas <- list()

for(i in 1:100){
  
  j <- i/100
  
  threshold <- hops_gdk3_ds1 > j
  
  num_cells <- cellStats(threshold, "sum")
  
  all_threshold_areas[i] <- num_cells
  
}

gdk3ds1 <- as.data.frame(all_threshold_areas)
gdk3ds1 <- rbind(gdk3ds1, rep(1:100))
gdk3ds1 <- t(gdk3ds1)
colnames(gdk3ds1) <- c("area", "threshold_percent")


#plotting

par(mar=c(5,4,4,4))

par(mfrow=c(4,2))
plot(hops_gdk0_ds10, main = c("GDK0 DS 10"))
plot(hops_gdk0_ds1, main = c("GDK0 DS 1"))

plot(hops_gdk1_ds10, main = c("GDK1 DS 10"))
plot(hops_gdk1_ds1, main = c("GDK1 DS 1"))

plot(hops_gdk3_ds10, main = c("GDK3 DS 10"))
plot(hops_gdk3_ds1, main = c("GDK3 DS 1"))



#gdk3
plot(gdk3ds10[,2] ~ gdk3ds10[,1], ylab = "Threshold", xlab = "Total pixels", main = c("Threshold Area Plot"), type="l", col = "black", lty = 1, xlim = c(0, 3e+07))
lines(gdk3ds1[,2] ~ gdk3ds1[,1], type="l", col = "black", lty = 5)

#gdk1
lines(gdk1ds10[,2]~ gdk1ds10[,1], type = "l", col='blue', lty = 1)
lines(gdk1ds1[,2] ~ gdk1ds1[,1], type="l", col = "blue", lty = 5)

#gdk0
lines(gdk0ds10[,2]~ gdk0ds10[,1], type = "l", col='red', lty = 1)
lines(gdk0ds1[,2] ~ gdk0ds1[,1], type="l", col = "red", lty = 5)


legend("topright", legend=c("GDK3 DS10", "GDK3 DS1", "GDK1 DS10", "GDK1 DS1", "GDK0 DS10", "GDK0 DS1"), col = c("black", "black", "blue", "blue", "red", "red"), lty = c(1,5,1,5,1,5))

#hops

bias_list <- get(paste(species.names[7], "_bias_list", sep="")) #get the species_bias_list object 

#prepare species occurrence and background data
d <- prepareData(species.file.paths[7], downsample_value = downsample.values[1], gdk_value = gdk.values[1], bias_file = bias_list[1])

plot(mean.temp.bio1, col="gray", legend=FALSE, main="hops Occurrences")
points(d[[1]], pch = 19)





####################################################################################################################################################################################################################

####################################################################################################################################################################################################################

####################################################################################################################################################################################################################


# Script for predicting all maxent .rds files 

# June 30 2019

#toadflax

setwd("F:/ENMEval Models")

base <- "F:/ENMEval Models/"
folders <- list.files()
fp <- paste(base, folders[3], sep="")
rda_files <- list.files(fp, full.names = TRUE)

for(i in 1:length(rda_files)){
  
  rda_file <- load(rda_files[i])
  maxent.prediction <- rmaxent::project(best.model, bioclim.jh.3)
  name <- paste(substr(rda_files[i], 1, nchar(rda_files[i])-4), ".tif", sep="")
  writeRaster(maxent.prediction$prediction_cloglog, filename = name, overwrite = TRUE)
  gc()
  
}

























































