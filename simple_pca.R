#script for downloading, processing, 
#PRISM climate data for common tansy 
#climate change common garden experiment

library(devtools)
library(prism)
library(ggplot2)


#set output folder
prism_set_dl_dir("D:/Common Tansy Climate Analysis/temp monthly norms")

get_prism_normals(type="tmean", resolution = "800m", mon = 1:12, keepZip=TRUE)

prism_set_dl_dir("D:/Common Tansy Climate Analysis/precip monthly norms")

get_prism_normals(type="ppt", resolution = "800m", mon = 1:12, keepZip=TRUE)


#for each lat, long site, extract 


tansy_climate  <- read.csv("C:/users/thomas/desktop/Tansy Seed Collections Master Climate.csv")

ncol(tansy_climate)

s <- cov(tansy_climate[,2:25])

sum(diag(s))

s.eigen <- eigen(s)
s.eigen

#first two PC account for 75% of variance
for(s in s.eigen$values){ print(s/sum(s.eigen$values))}

#scree graph of eigenvalue size vs eigenvalue number
plot(s.eigen$values, xlab = 'Eigenvalue Number', ylab = 'Eigenvalue Size', main = 'Scree Graph')
lines(s.eigen$values)

tansy.pca <- prcomp(tansy_climate[,5:ncol(tansy_climate)])

library(ggfortify)

pca.plot.scaled <- autoplot(tansy.pca, data = tansy_climate, colour = 'Record', label = TRUE)
pca.plot.scaled

kmeans <- autoplot(kmeans(tansy_climate[2:25], 6), data = tansy_climate, label=TRUE, label.size=4)
kmeans

library(cluster)
autoplot(fanny(tansy_climate[2:25], 6), frame=TRUE)


precip <- read.csv("c:/users/thomas/desktop/mn_climate/precip_climate.csv")

precip <- precip[,2:4]

precip_wide <- reshape(precip, idvar = "SrcID_Feat", v.names="Value", direction = "wide")









spurge_pca <- read.csv("c:/users/thomas/desktop/spurge_pca_aoi2_2020.csv")

s <- cov(spurge_pca[,3:10])

sum(diag(s))

s.eigen <- eigen(s)

for(s in s.eigen$values){ print(s/sum(s.eigen$values))}

#scree graph of eigenvalue size vs eigenvalue number
plot(s.eigen$values, xlab = 'Eigenvalue Number', ylab = 'Eigenvalue Size', main = 'Scree Graph')
lines(s.eigen$values)

spurge.pca <- prcomp(spurge_pca[,2:9])

library(ggfortify)

pca.plot.scaled <- autoplot(spurge.pca, data = spurge_pca, colour = 'Class', label = FALSE)
pca.plot.scaled



#satellite image t-SNE

#install.packages("tsne")
library(Rtsne)
library(tsne)
library(tidyverse)
library(rio)
library(RColorBrewer)
library(knitr)
library(plotly)

spurge_pca <- read.csv("c:/users/thomas/desktop/spurge_pca_aoi2_2020.csv")

Labels <- spurge_pca$Class

colors = rainbow(length(unique(Labels)))
names(colors) = unique(Labels)
## Executing the algorithm on curated data
tsne <- Rtsne(spurge_pca[,2:9], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500,check_duplicates=FALSE)


colourCount = length(unique(spurge_pca$Class))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
df=data_frame(Y1=tsne$Y[,1],Y2=tsne$Y[,2],label=Labels)
plot(df$Y1,df$Y2,col=colors[Labels], main="tsne")

ggplot(df,aes(Y1,Y2,colour=label))+geom_point()





















