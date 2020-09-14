###########################################################
# Supplementary Bivariate LiSA Figure
# Created by: Hari Iyer 
# Date: 4/21/2020
###########################################################


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

tz_shp@data$clust_HH_2 <-tz_shp@data$go_fdr_HH
tz_shp@data$clust_LL_2 <-tz_shp@data$go_fdr_LL
tz_shp@data$clust_LH_2 <-tz_shp@data$go_fdr_LH
tz_shp@data$clust_HL_2 <-tz_shp@data$go_fdr_HL

## functions

getCluster<-function(indat){
  indat$cluster<-0
  indat$cluster[indat$clust_HH_2==1]<-1
  indat$cluster[indat$clust_LL_2==1]<-2
  indat$cluster[indat$clust_LH_2==1]<-3
  indat$cluster[indat$clust_HL_2==1]<-4
  
  indat$cluster.f<-factor(indat$cluster, levels = c("0", "1", "2","3","4"), labels=c("Not significant","High/Long","Low/Short","Low/Long","High/Short"))
  ## 0 = NS, 1 = HH, 2 = LL, 3=LH, 4=HL
  return(indat)
  
}


## make maps

library(RColorBrewer)
#mycols<-c("gray57","firebrick4", "blue4","skyblue3","chartreuse4")
#mycols_k<-c("gray57","firebrick4","skyblue3","chartreuse4")
#mycols<-c("gray78","firebrick4", "darkolivegreen1","indianred1","darkgreen")
mycols<-c("gray78","#440154FF","#55C667FF","#39568CFF","#FDE725FF")

ken_shp@data<-getCluster(ken_shp@data)
p1 <- tm_shape(ken_shp) + tm_polygons(col = "cluster.f", palette=mycols, style="cat", title="Clusters", border.col = "gray18") + tm_layout("A", legend.position=c("left","bottom"))+tm_compass(position=c("right","top")) + tm_scale_bar(width=0.2) 

tz_shp@data<-getCluster(tz_shp@data)
p2<-tm_shape(tz_shp) + tm_polygons(col = "cluster.f", palette=mycols, style="cat", title="Clusters", border.col = "gray18") + tm_layout("B", legend.position=c("left","bottom")) +tm_compass(position=c("right","top")) + tm_scale_bar(width=0.15)


rwa_shp@data<-getCluster(rwa_shp@data)
p3<-tm_shape(rwa_shp) + tm_polygons(col = "cluster.f", palette=mycols, style="cat", title="Clusters", border.col = "gray18") + tm_layout("C", legend.position=c("left","top")) +tm_compass(position=c("right","top")) + tm_scale_bar(width=0.2)


mwi_shp@data<-getCluster(mwi_shp@data)
p4<-tm_shape(mwi_shp) + tm_polygons(col = "cluster.f", palette=mycols, style="cat", title="Clusters", border.col = "gray18") + tm_layout("D", legend.position=c("left","bottom")) +tm_compass(position=c("right","top")) + tm_scale_bar(width=0.2)


tmap_arrange(p1,p2,p3,p4, nrow=2)