###libraries####
library(ggplot2)
library(cowplot)
library(broom)
# library(lattice)
library(car)
library(reshape2)
# library(reshape)
# library(TSA)
library(plyr)
library(dplyr)
library(tidyverse)
# library(tidyr)
# library(rmarkdown)
library(visreg)
# library(MASS)
library(modEvA)
# library(BiodiversityR)
# library(gridExtra)
# library(AICcmodavg)
# library(nlme)
# library(mgcv)
# library(lme4)
# library(lsmeans)
library(multcomp)

# library(VTrack)
# library(igraph)
# library(MARSS)
library(splitstackshape) ###package to use cSplit and do text to column like excel
library(chron)
library(lubridate)
# library(rgdal)#package for geospatial data
# library(RInSp)#package for intraspecific niche variation
# library(boot)#package for boostrapping operations

############################################################
#Additional Packages
library(simmr)
library(SIBER)
library(viridis)
library(MixSIAR)

######################Data###########################33



###
#Pinfish mixing models results done by zone
###
pin.mix<-read.csv("BBpinfish.csv")

mix.z1 = subset(pin.mix, Zone == "Zone 1")
mix.z2 = subset(pin.mix, Zone == "Zone 2")

###
#Source data by zone
###
source_z1 = read.csv("sz1.csv")
source_z2 = read.csv("sz2.csv")

#################################################################



#################Data Processing/Preparation#################################

# calculate z scores----
###

#Zone 1
mix.z1$zBenthic.Algae = (mix.z1$Benthic.algae - mean(mix.z1$Benthic.algae))/sd(mix.z1$Benthic.algae)
mix.z1$zEpiphytes = (mix.z1$Epiphytes - mean(mix.z1$Epiphytes))/sd(mix.z1$Epiphytes)
mix.z1$zSeagrass = (mix.z1$Seagrass - mean(mix.z1$Seagrass))/sd(mix.z1$Seagrass)
mix.z1$zDrift = (mix.z1$Drift - mean(mix.z1$Drift))/sd(mix.z1$Drift)
mix.z1$zTL = (mix.z1$TL - mean(mix.z1$TL))/sd(mix.z1$TL)

#Zone 2
mix.z2$zBenthic.Algae = (mix.z2$Benthic.algae - mean(mix.z2$Benthic.algae))/sd(mix.z2$Benthic.algae)
mix.z2$zEpiphytes = (mix.z2$Epiphytes - mean(mix.z2$Epiphytes))/sd(mix.z2$Epiphytes)
mix.z2$zSeagrass = (mix.z2$Seagrass - mean(mix.z2$Seagrass))/sd(mix.z2$Seagrass)
mix.z2$zDrift = (mix.z2$Drift - mean(mix.z2$Drift))/sd(mix.z2$Drift)
mix.z2$zTL = (mix.z2$TL - mean(mix.z2$TL))/sd(mix.z2$TL)

#########################################################################################################################



###############Hypervolume estimation######################################################################

# generate hypervolumes----
library(hypervolume)
#subset data 
#use this to subset per site to calculate hypervolume site specific

#Zone 1
pin.s48 <- filter(mix.z1, Site == 48)%>%dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)
pin.s94 <- filter(mix.z1, Site == 94)%>%dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)
pin.s80 <- filter(mix.z1, Site == 80)%>%dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)
pin.s64 <- filter(mix.z1, Site == 64)%>%dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)
pin.s74 <- filter(mix.z1, Site == 74)%>%dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)
pin.s78 <- filter(mix.z1, Site == 78)%>%dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)

#Zone 2
pin.s12 <- filter(mix.z2, Site == 12)%>%dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)
pin.s23 <- filter(mix.z2, Site == 23)%>%dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)
pin.s07 <- filter(mix.z2, Site == 07)%>%dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)
pin.s18 <- filter(mix.z2, Site == 18)%>%dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)
pin.s13 <- filter(mix.z2, Site == 13)%>%dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)
pin.s02 <- filter(mix.z2, Site == 02)%>%dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)

# function to bootstrap hypervolumes to generate confidence intervals----
# makes "count" number of hvs by resampling number of points 
# equal to a fractional value or number of samples (percent)
# from dataset "values" and stores volumes

#re-sampling of hypervolumes using a subset of data
makerand<-function(values,count,percent){
  df_tot = data.frame(Volume = integer(0))
  number = ceiling(percent*nrow(values))
  if(percent>1){
    number = percent
  }
  for (i in 1:count){
    require(dplyr)
    rand = sample_n(values, number, replace = TRUE)
    require(hypervolume)
    randh = hypervolume_gaussian(rand,
                                 samples.per.point = ceiling((10^(3 + sqrt(ncol(rand))))/nrow(rand)),
                                 kde.bandwidth = estimate_bandwidth(rand), 
                                 sd.count = 3, 
                                 quantile.requested = 0.95, 
                                 quantile.requested.type = "probability", 
                                 chunk.size = 1000, 
                                 verbose = TRUE)
    volume=as.numeric(get_volume(randh))
    df=data.frame(Volume=volume)
    df_tot = rbind(df_tot, df)
    print(df_tot)
  }
  return(df_tot)
}

# function to calculate CI from bootstapped hypervolumes----
# Returns a CI on the hvs (WORKS IF DATASET IS 1000 or 100)
# hvconfint<-function(d, interval, name){
#   low = ((50-(interval/2))/100)*(nrow(d))
#   high = ((50+(interval/2))/100)*(nrow(d))
#   dn = d[order(d$sorenson),]
#   lci = dn$sorenson[low]
#   uci = dn$sorenson[high]
#   df=data.frame(Name=as.factor(name),Lower=lci,Upper=uci)
# }

#Change this to fix error and produce 95% based on the boostrapped results
hvconfint<-function(d, interval, name){
  low = ((50-(interval/2))/100)*(nrow(d))
  high = ((50+(interval/2))/100)*(nrow(d))
  dn = d[order(d$Volume),]
  lci = dn$Volume[low]
  uci = dn$Volume[high]
  df=data.frame(Name=as.factor(name),Lower=lci,Upper=uci)
}

###########################################################
##Pinfish Hypervolumes Zone 1----

# Pinfish Site 48 make hypervolume
pin.s48h = hypervolume_gaussian(pin.s48, name = 'Pinfish Site 48', 
                               samples.per.point = ceiling((10^(3 + sqrt(ncol(pin.s48))))/nrow(pin.s48)),
                               kde.bandwidth = estimate_bandwidth(pin.s48), 
                               sd.count = 3, 
                               quantile.requested = 0.95, 
                               quantile.requested.type = "probability", 
                               chunk.size = 1000, 
                               verbose = TRUE)
# get volume
get_volume(pin.s48h)

# # generate 1000 bootstrapped hvs with 9/10 of data - I changed this to 9/10 due the sample size (N = 10) per site
# pin.s48fish = makerand(pin.s48,1000, (9/10))
# pin.s48fish$cat = "Pinfish Site 48"
# 
# # calculate 95% confidence intervals from bootstraped hvs
# pin.s48conf = hvconfint(pin.s48fish, 95, "Pinfish Site 48")

# plot pinfish site 48 hypervolume

tiff("hypervol.pin.s48.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(pin.s48h, pairplot = T, colors=c("firebrick"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=F, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()



# Pinfish Site 94 make hypervolume
pin.s94h = hypervolume_gaussian(pin.s94, name = 'Pinfish Site 94', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(pin.s94))))/nrow(pin.s94)),
                                kde.bandwidth = estimate_bandwidth(pin.s94), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)
# get volume
get_volume(pin.s94h)

# # generate 1000 bootstrapped hvs with 9/10 of data
# pin.s94fish = makerand(pin.s94,1000, (9/10))
# pin.s94fish$cat = "Pinfish Site 94"
# 
# # calculate 95% confidence intervals from bootstraped hvs
# pin.s94conf = hvconfint(pin.s94fish, 95, "Pinfish Site 94")

# plot pinfish site 94 hypervolume

tiff("hypervol.pin.s94.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(pin.s94h, pairplot = T, colors=c("firebrick"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=F, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()




# Pinfish Site 80 make hypervolume
pin.s80h = hypervolume_gaussian(pin.s80, name = 'Pinfish Site 80', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(pin.s80))))/nrow(pin.s80)),
                                kde.bandwidth = estimate_bandwidth(pin.s80), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)
# get volume
get_volume(pin.s80h)

# # generate 1000 bootstrapped hvs with 9/10 of data
# pin.s80fish = makerand(pin.s80,1000, (9/10))
# pin.s80fish$cat = "Pinfish Site 80"
# 
# # calculate 95% confidence intervals from bootstraped hvs
# pin.s80conf = hvconfint(pin.s80fish, 95, "Pinfish Site 80")

# plot pinfish site 80 hypervolume

tiff("hypervol.pin.s80.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(pin.s80h, pairplot = T, colors=c("firebrick"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=F, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()



# Pinfish Site 64 make hypervolume
pin.s64h = hypervolume_gaussian(pin.s64, name = 'Pinfish Site 64', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(pin.s64))))/nrow(pin.s64)),
                                kde.bandwidth = estimate_bandwidth(pin.s64), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)
# get volume
get_volume(pin.s64h)

# # generate 1000 bootstrapped hvs with 9/10 of data
# pin.s64fish = makerand(pin.s64,1000, (9/10))
# pin.s64fish$cat = "Pinfish Site 64"
# 
# # calculate 95% confidence intervals from bootstraped hvs
# pin.s64conf = hvconfint(pin.s64fish, 95, "Pinfish Site 64")

# plot pinfish site 64 hypervolume

tiff("hypervol.pin.s64.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(pin.s64h, pairplot = T, colors=c("firebrick"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=F, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()


# Pinfish Site 74 make hypervolume
pin.s74h = hypervolume_gaussian(pin.s74, name = 'Pinfish Site 74', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(pin.s74))))/nrow(pin.s74)),
                                kde.bandwidth = estimate_bandwidth(pin.s74), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)
# get volume
get_volume(pin.s74h)
# 
# # generate 1000 bootstrapped hvs with 9/10 of data
# pin.s74fish = makerand(pin.s74,1000, (9/10))
# pin.s74fish$cat = "Pinfish Site 74"
# 
# # calculate 95% confidence intervals from bootstraped hvs
# pin.s74conf = hvconfint(pin.s74fish, 95, "Pinfish Site 74")

# plot pinfish site 74 hypervolume

tiff("hypervol.pin.s74.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(pin.s74h, pairplot = T, colors=c("firebrick"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=F, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()



# Pinfish Site 78 make hypervolume
pin.s78h = hypervolume_gaussian(pin.s78, name = 'Pinfish Site 78', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(pin.s78))))/nrow(pin.s78)),
                                kde.bandwidth = estimate_bandwidth(pin.s78), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)
# get volume
get_volume(pin.s78h)

# # generate 1000 bootstrapped hvs with 9/10 of data
# pin.s78fish = makerand(pin.s78,1000, (9/10))
# pin.s78fish$cat = "Pinfish Site 78"
# 
# # calculate 95% confidence intervals from bootstraped hvs
# pin.s78conf = hvconfint(pin.s78fish, 95, "Pinfish Site 78")

# plot pinfish site 78 hypervolume

tiff("hypervol.pin.s78.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(pin.s78h, pairplot = T, colors=c("firebrick"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=F, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()


########################################################################################################

####
#Adding volume for all sites

pin.vol.z1 <- as.data.frame(cbind(get_volume(pin.s48h), get_volume(pin.s94h), get_volume(pin.s80h),
                    get_volume(pin.s64h), get_volume(pin.s74h), get_volume(pin.s78h)))
colnames(pin.vol.z1)<-c(48, 94, 80, 64, 74, 78)
rownames(pin.vol.z1)<-c("Volume")
pin.vol.z1 <- gather(pin.vol.z1, "Site_Num", "Volume")
pin.vol.z1$Zone <- c("Zone 1") 

###
#Bootstrapped hvs Zone 1 -----


# generate 1000 bootstrapped hvs with 9/10 of data - I changed this to 9/10 due the sample size (N = 10) per site
pin.s48fish = makerand(pin.s48,1000, (9/10))
pin.s48fish$cat = "Pinfish Site 48"

# calculate 95% confidence intervals from bootstraped hvs
pin.s48conf = hvconfint(pin.s48fish, 95, "Pinfish Site 48")

# generate 1000 bootstrapped hvs with 9/10 of data
pin.s80fish = makerand(pin.s80,1000, (9/10))
pin.s80fish$cat = "Pinfish Site 80"

# calculate 95% confidence intervals from bootstraped hvs
pin.s80conf = hvconfint(pin.s80fish, 95, "Pinfish Site 80")

# generate 1000 bootstrapped hvs with 9/10 of data
pin.s94fish = makerand(pin.s94,1000, (9/10))
pin.s94fish$cat = "Pinfish Site 94"

# calculate 95% confidence intervals from bootstraped hvs
pin.s94conf = hvconfint(pin.s94fish, 95, "Pinfish Site 94")

# generate 1000 bootstrapped hvs with 9/10 of data
pin.s64fish = makerand(pin.s64,1000, (9/10))
pin.s64fish$cat = "Pinfish Site 64"

# calculate 95% confidence intervals from bootstraped hvs
pin.s64conf = hvconfint(pin.s64fish, 95, "Pinfish Site 64")

# generate 1000 bootstrapped hvs with 9/10 of data
pin.s74fish = makerand(pin.s74,1000, (9/10))
pin.s74fish$cat = "Pinfish Site 74"

# calculate 95% confidence intervals from bootstraped hvs
pin.s74conf = hvconfint(pin.s74fish, 95, "Pinfish Site 74")

# generate 1000 bootstrapped hvs with 9/10 of data
pin.s78fish = makerand(pin.s78,1000, (9/10))
pin.s78fish$cat = "Pinfish Site 78"

# calculate 95% confidence intervals from bootstraped hvs
pin.s78conf = hvconfint(pin.s78fish, 95, "Pinfish Site 78")




##########################################################################################################
##Pinfish Hypervolumes Zone 2----

# Pinfish Site 12 make hypervolume
pin.s12h = hypervolume_gaussian(pin.s12, name = 'Pinfish Site 12', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(pin.s12))))/nrow(pin.s12)),
                                kde.bandwidth = estimate_bandwidth(pin.s12), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)
# get volume
get_volume(pin.s12h)

# # generate 1000 bootstrapped hvs with 9/10 of data
# pin.s12fish = makerand(pin.s12,1000, (9/10))
# pin.s12fish$cat = "Pinfish Site 12"
# 
# # calculate 95% confidence intervals from bootstraped hvs
# pin.s12conf = hvconfint(pin.s12fish, 95, "Pinfish Site 12")

# plot pinfish site 12 hypervolume

tiff("hypervol.pin.s12.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(pin.s12h, pairplot = T, colors=c("firebrick"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=F, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()



# Pinfish Site 23 make hypervolume
pin.s23h = hypervolume_gaussian(pin.s23, name = 'Pinfish Site 23', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(pin.s23))))/nrow(pin.s23)),
                                kde.bandwidth = estimate_bandwidth(pin.s23), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)
# get volume
get_volume(pin.s23h)

# generate 1000 bootstrapped hvs with 9/10 of data
# pin.s23fish = makerand(pin.s23,1000, (9/10))
# pin.s23fish$cat = "Pinfish Site 23"
# 
# # calculate 95% confidence intervals from bootstraped hvs
# pin.s23conf = hvconfint(pin.s23fish, 95, "Pinfish Site 23")

# plot pinfish site 23 hypervolume

tiff("hypervol.pin.s23.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(pin.s23h, pairplot = T, colors=c("firebrick"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=F, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()




# Pinfish Site 07 make hypervolume
pin.s07h = hypervolume_gaussian(pin.s07, name = 'Pinfish Site 07', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(pin.s07))))/nrow(pin.s07)),
                                kde.bandwidth = estimate_bandwidth(pin.s07), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)
# get volume
get_volume(pin.s07h)

# # generate 1000 bootstrapped hvs with 9/10 of data
# pin.s07fish = makerand(pin.s07,1000, (9/10))
# pin.s07fish$cat = "Pinfish Site 07"
# 
# # calculate 95% confidence intervals from bootstraped hvs
# pin.s07conf = hvconfint(pin.s07fish, 95, "Pinfish Site 07")

# plot pinfish site 07 hypervolume

tiff("hypervol.pin.s07.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(pin.s07h, pairplot = T, colors=c("firebrick"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=F, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()



# Pinfish Site 18 make hypervolume
pin.s18h = hypervolume_gaussian(pin.s18, name = 'Pinfish Site 18', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(pin.s18))))/nrow(pin.s18)),
                                kde.bandwidth = estimate_bandwidth(pin.s18), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)
# get volume
get_volume(pin.s18h)

# # generate 1000 bootstrapped hvs with 9/10 of data
# pin.s18fish = makerand(pin.s18,1000, (9/10))
# pin.s18fish$cat = "Pinfish Site 18"
# 
# # calculate 95% confidence intervals from bootstraped hvs
# pin.s18conf = hvconfint(pin.s18fish, 95, "Pinfish Site 18")

# plot pinfish site 18 hypervolume

tiff("hypervol.pin.s18.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(pin.s18h, pairplot = T, colors=c("firebrick"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=F, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()


# Pinfish Site 13 make hypervolume
pin.s13h = hypervolume_gaussian(pin.s13, name = 'Pinfish Site 13', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(pin.s13))))/nrow(pin.s13)),
                                kde.bandwidth = estimate_bandwidth(pin.s13), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)
# get volume
get_volume(pin.s13h)

# # generate 1000 bootstrapped hvs with 9/10 of data
# pin.s13fish = makerand(pin.s13,1000, (9/10))
# pin.s13fish$cat = "Pinfish Site 13"
# 
# # calculate 95% confidence intervals from bootstraped hvs
# pin.s13conf = hvconfint(pin.s13fish, 95, "Pinfish Site 13")

# plot pinfish site 13 hypervolume

tiff("hypervol.pin.s13.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(pin.s13h, pairplot = T, colors=c("firebrick"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=F, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()



# Pinfish Site 02 make hypervolume
pin.s02h = hypervolume_gaussian(pin.s02, name = 'Pinfish Site 02', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(pin.s02))))/nrow(pin.s02)),
                                kde.bandwidth = estimate_bandwidth(pin.s02), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)
# get volume
get_volume(pin.s02h)

# # generate 1000 bootstrapped hvs with 9/10 of data
# pin.s02fish = makerand(pin.s02,1000, (9/10))
# pin.s02fish$cat = "Pinfish Site 02"
# 
# # calculate 95% confidence intervals from bootstraped hvs
# pin.s02conf = hvconfint(pin.s02fish, 95, "Pinfish Site 02")

# plot pinfish site 02 hypervolume

tiff("hypervol.pin.s02.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(pin.s02h, pairplot = T, colors=c("firebrick"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=F, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()


########################################################################################################

####
#Adding volume for all sites

pin.vol.z2 <- as.data.frame(cbind(get_volume(pin.s12h), get_volume(pin.s23h), get_volume(pin.s07h),
                                  get_volume(pin.s18h), get_volume(pin.s13h), get_volume(pin.s02h)))
colnames(pin.vol.z2)<-c(12, 23, 07, 18, 13, 02)
rownames(pin.vol.z2)<-c("Volume")
pin.vol.z2 <- gather(pin.vol.z2, "Site_Num", "Volume")
pin.vol.z2$Zone <- c("Zone 2") 


write.csv(pin.vol.z1, "pin.vol.z1.csv")
write.csv(pin.vol.z2, "pin.vol.z2.csv")

pin.vol.all <-bind_rows(pin.vol.z1, pin.vol.z2)
write.csv(pin.vol.all, "pin.vol.all.csv")

###
#Bootstrapped hvs Zone 2 -----

# generate 1000 bootstrapped hvs with 9/10 of data - I changed this to 9/10 due the sample size (N = 10) per site
pin.s12fish = makerand(pin.s12,1000, (9/10))
pin.s12fish$cat = "Pinfish Site 12"

# calculate 95% confidence intervals from bootstraped hvs
pin.s12conf = hvconfint(pin.s12fish, 95, "Pinfish Site 12")

# generate 1000 bootstrapped hvs with 9/10 of data
pin.s23fish = makerand(pin.s23,1000, (9/10))
pin.s23fish$cat = "Pinfish Site 23"

# calculate 95% confidence intervals from bootstraped hvs
pin.s23conf = hvconfint(pin.s23fish, 95, "Pinfish Site 23")

# generate 1000 bootstrapped hvs with 9/10 of data
pin.s07fish = makerand(pin.s07,1000, (9/10))
pin.s07fish$cat = "Pinfish Site 07"

# calculate 95% confidence intervals from bootstraped hvs
pin.s07conf = hvconfint(pin.s07fish, 95, "Pinfish Site 07")

# generate 1000 bootstrapped hvs with 9/10 of data
pin.s18fish = makerand(pin.s18,1000, (9/10))
pin.s18fish$cat = "Pinfish Site 18"

# calculate 95% confidence intervals from bootstraped hvs
pin.s18conf = hvconfint(pin.s18fish, 95, "Pinfish Site 18")

# generate 1000 bootstrapped hvs with 9/10 of data
pin.s13fish = makerand(pin.s13,1000, (9/10))
pin.s13fish$cat = "Pinfish Site 13"

# calculate 95% confidence intervals from bootstraped hvs
pin.s13conf = hvconfint(pin.s13fish, 95, "Pinfish Site 13")

# generate 1000 bootstrapped hvs with 9/10 of data
pin.s02fish = makerand(pin.s02,1000, (9/10))
pin.s02fish$cat = "Pinfish Site 02"

# calculate 95% confidence intervals from bootstraped hvs
pin.s02conf = hvconfint(pin.s02fish, 95, "Pinfish Site 02")



##########################################################################


##################Overlap analysis#####################################

###
#Prepare data to create general hypervolume - per zone, and per zone|seascape---------
###

zone.1 <- filter(mix.z1, Zone %in% c("Zone 1"))%>%dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)
zone.2 <- filter(mix.z2, Zone %in% c("Zone 2"))%>%dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)

zone.1.c <- filter(mix.z1, Zone %in% c("Zone 1") & Site.y %in% c(94, 80, 48))%>%
  dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)

zone.1.f <- filter(mix.z1, Zone %in% c("Zone 1") & Site.y %in% c(64, 74, 78))%>%
  dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)

zone.2.c <- filter(mix.z2, Zone %in% c("Zone 2") & Site.y %in% c(23, 12, 7))%>%
  dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)

zone.2.f <- filter(mix.z2, Zone %in% c("Zone 2") & Site.y %in% c(18, 13, 2))%>%
  dplyr::select(.,zBenthic.Algae, zEpiphytes, zSeagrass, zDrift, zTL)

###

#Function from Ryan/Lesser
###
# function to generate bootstapped overlaps between 2 hypervolumes----
# overlap stats - make a hv of fraction of random samples of set 1, 
# make a hv of fraction of random samples of set 2, get and store overlap, repeat count times
randoverlap<-function(data1,data2,count,percent){
  require(hypervolume)
  df_tot = data.frame(sorenson = integer(0), unique1 = integer(0), unique2 = integer(0))
  number1 = ceiling(percent*nrow(data1))
  if(percent>1){
    number1 = percent
  }
  number2 = ceiling(percent*nrow(data2))
  if(percent>1){
    number2 = percent
  }
  for (i in 1:count){
    require(dplyr)
    rand1 = sample_n(data1, number1, replace = TRUE)
    rand2 = sample_n(data2, number2, replace = TRUE)
    require(hypervolume)
    randh1 = hypervolume_gaussian(rand1,
                                  samples.per.point = ceiling((10^(3 + sqrt(ncol(rand1))))/nrow(rand1)),
                                  kde.bandwidth = estimate_bandwidth(rand1), 
                                  sd.count = 3, 
                                  quantile.requested = 0.95, 
                                  quantile.requested.type = "probability", 
                                  chunk.size = 1000, 
                                  verbose = TRUE)
    randh2 = hypervolume_gaussian(rand2,
                                  samples.per.point = ceiling((10^(3 + sqrt(ncol(rand2))))/nrow(rand2)),
                                  kde.bandwidth = estimate_bandwidth(rand2), 
                                  sd.count = 3, 
                                  quantile.requested = 0.95, 
                                  quantile.requested.type = "probability", 
                                  chunk.size = 1000, 
                                  verbose = TRUE)
    set<-hypervolume_set(randh1,randh2,check.memory = FALSE)
    ov<-data.frame(hypervolume_overlap_statistics(set))
    soi<-ov[2,]
    frac1<-ov[3,]
    frac2<-ov[4,]
    df=data.frame(sorenson=soi,unique1 = frac1, unique2 = frac2)
    df_tot = rbind(df_tot, df)
    print(df_tot)
    cat('done with ', i, '/',count,'\n')
  }
  return(df_tot)
}


## calculate overlap between Zones----

###
#calc volume within each zone
###

zone.1.h = hypervolume_gaussian(zone.1, name = 'Pinfish Zone 1', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(zone.1))))/nrow(zone.1)),
                                kde.bandwidth = estimate_bandwidth(zone.1), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)

zone.2.h = hypervolume_gaussian(zone.2, name = 'Pinfish Zone 2', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(zone.2))))/nrow(zone.2)),
                                kde.bandwidth = estimate_bandwidth(zone.2), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)

###
#Generate hypervolume figure Zone Overlap
###
hZoneAll <- hypervolume_join(zone.1.h, zone.2.h)

tiff("hypervol.ZonesOverlap.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(hZoneAll, pairplot = T, colors=c("firebrick", "#56B4E9"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=T, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()



###
#Overlap estimate and boostrap dist
###

zone.comp = hypervolume_set(zone.1.h, zone.2.h, check.memory = FALSE)

# generate overlap
hypervolume_overlap_statistics(zone.comp)

# generate bootstap overlaps
zone.comp.boot = randoverlap(zone.1, zone.2, 100, 2/3)
zone.comp.boot$cat = as.factor("Zones")
write_csv(zone.comp.boot, 'Zonesbootoverlaps.csv')
# calculate 95% CI

zone.comp.bconf = hvconfint(zone.comp.boot, 95, "Zones")



## calculate overlap between seascapes----

###
#calc volume within each seascape|zone
###

cs.z1.h = hypervolume_gaussian(zone.1.c, name = 'Pinfish Continuous Z1', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(zone.1.c))))/nrow(zone.1.c)),
                                kde.bandwidth = estimate_bandwidth(zone.1.c), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)

cs.z2.h = hypervolume_gaussian(zone.2.c, name = 'Pinfish Continuous Z2', 
                                samples.per.point = ceiling((10^(3 + sqrt(ncol(zone.2.c))))/nrow(zone.2.c)),
                                kde.bandwidth = estimate_bandwidth(zone.2.c), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = TRUE)


fs.z1.h = hypervolume_gaussian(zone.1.f, name = 'Pinfish Fragmented Z1', 
                               samples.per.point = ceiling((10^(3 + sqrt(ncol(zone.1.f))))/nrow(zone.1.f)),
                               kde.bandwidth = estimate_bandwidth(zone.1.f), 
                               sd.count = 3, 
                               quantile.requested = 0.95, 
                               quantile.requested.type = "probability", 
                               chunk.size = 1000, 
                               verbose = TRUE)

fs.z2.h = hypervolume_gaussian(zone.2.f, name = 'Pinfish Fragmented Z2', 
                               samples.per.point = ceiling((10^(3 + sqrt(ncol(zone.2.f))))/nrow(zone.2.f)),
                               kde.bandwidth = estimate_bandwidth(zone.2.f), 
                               sd.count = 3, 
                               quantile.requested = 0.95, 
                               quantile.requested.type = "probability", 
                               chunk.size = 1000, 
                               verbose = TRUE)



###
#Generate hypervolume figure Zone Overlap
###
hSeascapeZ1 <- hypervolume_join(cs.z1.h, fs.z1.h)
hSeascapeZ2 <- hypervolume_join(cs.z2.h, fs.z2.h)


tiff("hypervol.SeascapeOverlap.Z1.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(hSeascapeZ1, pairplot = T, colors=c("green", "red"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=T, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()

tiff("hypervol.SeascapeOverlap.Z2.tiff", width = 12, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(hSeascapeZ2, pairplot = T, colors=c("green", "red"),
     names= c("Benthic Algae", "Epiphytes", "Seagrass", "Drift", "TL"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=T, limits=c(-5,5), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)

dev.off()

###
#Overlap estimate and boostrap dist
###

z1.seascape.comp = hypervolume_set(fs.z1.h, cs.z1.h, check.memory = FALSE)
z2.seascape.comp = hypervolume_set(fs.z2.h, cs.z2.h, check.memory = FALSE)


# generate overlap
hypervolume_overlap_statistics(z1.seascape.comp)
hypervolume_overlap_statistics(z2.seascape.comp)

# generate bootstap overlaps z1



# calculate 95% CI
zone.comp.bconf = hvconfint(zone.comp.boot, 95, "Zones")





# mid productivity
# combine mid pigfish and pinfish
mid = hypervolume_set(midpigh, midpinh, check.memory = FALSE)

# generate overlap
hypervolume_overlap_statistics(mid)

# generate bootstrap overlaps
midcomp = randoverlap(midpin,midpig,100,2/3)
midcomp$cat = as.factor("Mid")

# calculate 95% CI
ovmidconf = hvconfint(midcomp, 95, "Mid Productivity")


##high productivity
# combine high pigfish and pinfish
high = hypervolume_set(highpigh, highpinh, check.memory = FALSE)

# generate overlap
hypervolume_overlap_statistics(high)

# generate bootstrap overlaps
highcomp = randoverlap(highpin,highpig,100,2/3)
highcomp$cat = as.factor("High")

# calculate 95% CI
ovhighconf = hvconfint(highcomp, 95, "High Productivity")


# bind all overlaps together and confidence intervals together
overlaps = rbind(lowcomp,midcomp,highcomp)
ovconfints = rbind(ovlowconf,ovmidconf,ovhighconf)

# overlap statistics and plot
dunn.test(x=pinpigover$sorenson,g=pinpigover$cat, kw=TRUE, method="bonferroni")
plot(sorenson~cat, data=pinpigover)

