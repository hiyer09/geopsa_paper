#######################################################################################
# Program 04 - Regression coefficients and confidence intervals
# Goal: Get betas for log-transformed travel time, log-transformed population density
# for each country
# Output: Table
#######################################################################################


library(foreign)
library(sf)
library(raster)
library(spData)
library(readr)
library(rgdal)
library(grid)
library(ggplot2)
library(tmap)
library(MapGAM)
library(gtools)
library(gmodels)
library(spdep)
library(ggplot2)
library(spatialreg)

#######################################
## Step 1: Read in all datasets
#######################################

ken_shp<-readOGR(dsn = path.expand("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/data/Kenya"),
                 layer = "KEN_t2_exp", stringsAsFactors = FALSE)
#summary(ken_shp)
ken_shp$POPMEAN<-ken_shp$POPMEAN_1
ken_shp$LOGPOPMEAN<-ken_shp$LOGPOPME_1 
mwi_shp<-readOGR(dsn = path.expand("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/data/Malawi"),
                 layer = "mal_dist_join", stringsAsFactors = FALSE)
#summary(mwi_shp)

rwa_shp<-readOGR(dsn = path.expand("C:/Users/Hari/Documents/PIH/consulting/Rwanda"),
                 layer = "RWA_export_merge", stringsAsFactors = FALSE)
#summary(rwa_shp)

tz_shp<-readOGR(dsn = path.expand("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/TZ paper final/data"),
                layer = "tz_dist_join_v9", stringsAsFactors = FALSE)
#summary(tz_shp)


## spearman's correlation
# install.packages("DescTools")
library(DescTools)
c<-cor.test(ken_shp@data$POPMEAN, ken_shp@data$T2MEAN, method="spearman")

li<-list(ken_shp, mwi_shp, rwa_shp, tz_shp)

opcor<-as.matrix(ncol=3, nrow=4)
for (i in seq_along(li)){
  opcor[i]<-getCor(li[i])
} 


## Spearman's correlation
SpearmanRho(ken_shp@data$POPMEAN, ken_shp@data$T2MEAN, conf.level = 0.95)
SpearmanRho(mwi_shp@data$POPMEAN, mwi_shp@data$T2MEAN, conf.level = 0.95)
SpearmanRho(rwa_shp@data$POPMEAN, rwa_shp@data$T2MEAN, conf.level = 0.95)
SpearmanRho(tz_shp@data$POPMEAN, tz_shp@data$T2MEAN, conf.level = 0.95)

a<-cor.test(ken_shp@data$LOGPOPMEAN, ken_shp@data$LOGT2MEAN)
ci<-rbind(a$conf.int)

getCor<-function(indat){
  a<-cor.test(indat@data$LOGPOPMEAN, indat@data$LOGT2MEAN)
  pe<-a$estimate
  ci<-rbind(a$conf.int)
  f<-cbind(pe,ci)
  return(f)
}

getCor(rwa_shp)
getCor(ken_shp)


li<-list(ken_shp, mwi_shp, rwa_shp, tz_shp)

opcor<-as.data.frame(matrix(,ncol=3, nrow=4))
for (i in seq_along(li)){
  opcor[i,]<-getCor(li[[i]])
} 

## spearman

getSpRho<-function(indat){
  SpearmanRho(indat@data$POPMEAN, indat@data$T2MEAN, conf.level = 0.95)
}

opspcor<-as.data.frame(matrix(,ncol=3, nrow=4))
for (i in seq_along(li)){
  opspcor[i,]<-getSpRho(li[[i]])
} 

SpearmanRho(ken_shp@data$POPMEAN, ken_shp@data$T2MEAN, conf.level = 0.95)
SpearmanRho(mwi_shp@data$POPMEAN, mwi_shp@data$T2MEAN, conf.level = 0.95)
SpearmanRho(rwa_shp@data$POPMEAN, rwa_shp@data$T2MEAN, conf.level = 0.95)
SpearmanRho(tz_shp@data$POPMEAN, tz_shp@data$T2MEAN, conf.level = 0.95)






## functions

getRegSet<-function(indat){
  indat<-indat@data[c("LOGPOPMEAN","LOGT2MEAN")]
  return(indat)
}


## process data

kendat<-getRegSet(ken_shp)
mwidat<-getRegSet(mwi_shp)
rwadat<-getRegSet(rwa_shp)
tzdat<-getRegSet(tz_shp)


getRegOp<-function(indat){
  
  ## run regression
  reg<-lm(formula = LOGT2MEAN ~ LOGPOPMEAN, data=indat)
  s<-summary(reg)
  
  # get betas
  betas <-as.data.frame(s["coefficients"])
  betas<- betas$coefficients.Estimate
  
  # get CIs
  conf<-as.data.frame(confint(reg))[2,]
  
  
  ## for a 20% increase
  a<-as.data.frame(0.8^betas[2])
  b<-0.8^conf
  
  c<-cbind(a,b)
  
}


kenop<-getRegOp(kendat)
mwiop<-getRegOp(mwidat)
rwaop<-getRegOp(rwadat)
tzop<-getRegOp(tzdat)

## spatial lag

## shape
ken_shp<-readOGR(dsn = path.expand("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/data/Kenya"),
                 layer = "KEN_t2_exp", stringsAsFactors = FALSE)


## neighbors

kengal <- read.gal("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/data/Kenya/KEN_t2_exp.GAL",override.id=TRUE)
kenwts<-nb2listw(kengal)
summary(kenwts)

getNb<-function(indat){
  wts<-nb2listw(indat)
  return(wts)
}

mwigal<-read.gal("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/data/Malawi/mal_dist_join.gal",override.id=TRUE)
rwagal<-read.gal("C:/Users/Hari/Documents/PIH/consulting/Rwanda/RWA_export_merge.gal",override.id=TRUE)
tzgal<-read.gal("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/TZ paper final/data/tz_dist_join_v4.gal",override.id=TRUE)

mwiwts<-getNb(mwigal)
rwawts<-getNb(rwagal)
tzwts<-nb2listw(tzgal,zero.policy=TRUE)


## moran's I
moran.test(ken_shp$LOGPOPMEAN,kenwts,randomisation=FALSE, alternative="two.sided")


## spatial lag
lagken <- lagsarlm(LOGT2MEAN ~ LOGPOPMEAN, data=ken_shp,kenwts)
lagmwi <- lagsarlm(LOGT2MEAN ~ LOGPOPMEAN, data=mwi_shp,mwiwts)
lagken <- lagsarlm(LOGT2MEAN ~ LOGPOPMEAN, data=ken_shp,kenwts)
lagken <- lagsarlm(LOGT2MEAN ~ LOGPOPMEAN, data=ken_shp,kenwts)


getLagRegOp<-function(indat, inwts){
  
  ## run regression
  reg<-lagsarlm(LOGT2MEAN ~ LOGPOPMEAN, data=indat,inwts)
  s<-summary(reg)
  
  # get betas
  betas <-as.data.frame(s["coefficients"])
  betas<- betas$coefficients
  
  # get CIs
  conf<-as.data.frame(confint(reg))[3,]
  
  
  ## for a 20% increase
  a<-as.data.frame(0.8^betas[2])
  b<-0.8^conf
  
  c<-cbind(a,b)
  
}


getLagRegOpNA<-function(indat, inwts){
  
  ## run regression
  reg<-lagsarlm(LOGT2MEAN ~ LOGPOPMEAN, data=indat,inwts, zero.policy = TRUE)
  s<-summary(reg)
  
  # get betas
  betas <-as.data.frame(s["coefficients"])
  betas<- betas$coefficients
  
  # get CIs
  conf<-as.data.frame(confint(reg))[3,]
  
  
  ## for a 20% increase
  a<-as.data.frame(0.8^betas[2])
  b<-0.8^conf
  
  c<-cbind(a,b)
  
}

kenlop<-getLagRegOp(ken_shp, kenwts)
mwilop<-getLagRegOp(mwi_shp,mwiwts)
rwalop<-getLagRegOp(rwa_shp,rwawts)
tzlop<-getLagRegOpNA(tz_shp,tzwts)

rbind(kenop,kenlop)
rbind(mwiop,mwilop)
rbind(rwaop,rwalop)
rbind(tzop,tzlop)
