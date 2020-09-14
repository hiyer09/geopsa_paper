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
library(DescTools)

#######################################
## Step 1: Read in all datasets
#######################################

ken_shp<-readOGR(dsn = path.expand("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/Concept Paper/data/Data"),
                 layer = "ken_gojoin_v3", stringsAsFactors = FALSE)
summary(ken_shp)

mwi_shp<-readOGR(dsn = path.expand("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/Concept Paper/data/Data"),
                 layer = "MWI_go_join_v3", stringsAsFactors = FALSE)
summary(mwi_shp)

rwa_shp<-readOGR(dsn = path.expand("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/Concept Paper/data/Data"),
                 layer = "rwa_go_join_v3", stringsAsFactors = FALSE)
summary(rwa_shp)

tz_shp<-readOGR(dsn = path.expand("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/TZ paper final/data"),
                layer = "tz_dist_join_v9", stringsAsFactors = FALSE)
summary(tz_shp)

tz_shp@data$GOMEAN2<-tz_shp@data$GOMEAN
tz_shp@data$LOGGOMEAN2<-tz_shp@data$LOGGOMEAN

##############################################################################################
# ANALYSIS
##############################################################################################

### Pearson's Correlation


li<-list(ken_shp, mwi_shp, rwa_shp, tz_shp)

getCor<-function(indat){
  a<-cor.test(indat@data$LOGPOPMEAN, indat@data$LOGGOMEAN2)
  pe<-a$estimate
  ci<-rbind(a$conf.int)
  f<-cbind(pe,ci)
  return(f)
}



opcor<-as.data.frame(matrix(,ncol=3, nrow=4))
for (i in seq_along(li)){
  opcor[i,]<-getCor(li[[i]])
} 


## Spearman's Correlation


getSpRho<-function(indat){
  SpearmanRho(indat@data$POPMEAN, indat@data$GOMEAN2, conf.level = 0.95)
}

opspcor<-as.data.frame(matrix(,ncol=3, nrow=4))
for (i in seq_along(li)){
  opspcor[i,]<-getSpRho(li[[i]])
} 

## Regression


getRegSet<-function(indat){
  indat<-indat@data[c("LOGPOPMEAN","LOGGOMEAN2")]
  return(indat)
}


getRegOp<-function(indat){
  
  ## run regression
  reg<-lm(formula = LOGGOMEAN2 ~ LOGPOPMEAN, data=indat)
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

opreg <-as.data.frame(matrix(,ncol=3, nrow=4))
for(i in seq_along(li)){
  prep<-getRegSet(li[[i]])
  opreg[i,]<-getRegOp(prep)
}


## Spatial Lag Model

## neighbors

getNb<-function(indat){
  wts<-nb2listw(indat)
  return(wts)
}

kengal <- read.gal("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/Concept Paper/data/Data/ken_gojoin_v3.GAL",override.id=TRUE)
kenwts<-nb2listw(kengal)
summary(kenwts)



mwigal<-read.gal("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/Concept Paper/data/Data/MWI_go_join_v3.gal",override.id=TRUE)
rwagal<-read.gal("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/Concept Paper/data/Data/rwa_go_join_v3.gal",override.id=TRUE)
tzgal<-read.gal("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/Concept Paper/data/Data/tz_go_join_v4.gal",override.id=TRUE)

mwiwts<-getNb(mwigal)
rwawts<-getNb(rwagal)
tzwts<-nb2listw(tzgal,zero.policy=TRUE)


## moran's I
moran.test(ken_shp$LOGPOPMEAN,kenwts,randomisation=FALSE, alternative="two.sided")


## spatial lag
lagken <- lagsarlm(LOGGOMEAN2 ~ LOGPOPMEAN, data=ken_shp,kenwts)
lagmwi <- lagsarlm(LOGGOMEAN2 ~ LOGPOPMEAN, data=mwi_shp,mwiwts)
lagken <- lagsarlm(LOGGOMEAN2 ~ LOGPOPMEAN, data=ken_shp,kenwts)
lagken <- lagsarlm(LOGGOMEAN2 ~ LOGPOPMEAN, data=ken_shp,kenwts)


getLagRegOp<-function(indat, inwts){
  
  ## run regression
  reg<-lagsarlm(LOGGOMEAN2 ~ LOGPOPMEAN, data=indat,inwts)
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
  reg<-lagsarlm(LOGGOMEAN2 ~ LOGPOPMEAN, data=indat,inwts, zero.policy = TRUE)
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

rbind(kenlop,tzlop,rwalop,mwilop)
