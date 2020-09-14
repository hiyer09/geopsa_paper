#################################################################
# R Code for Bivariate LiSA Scatter - Concept Paper
# Objective: Evaluate urban-rural access to cancer research centers
# in SSA
# Goals: read in all GeoDA shapefiles
# Do the combined scatter plot, with factor variable for different
# cluster types


# Change the axis to reflect in-country - updated
# Update to tier 2
# Updated to incorporate FDR (B-H correction)

# Updated 9/3/2020 to respond to reviewers
# Changed x-axis label
# Changed colors
# Updated TZ dataset (w/500m resolution TT)
##################################################################


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
summary(ken_shp)
ken_shp = ken_shp[,!(names(ken_shp) %in% c("POPMEAN","LOGPOPMEAN"))]
ken_shp$POPMEAN<-ken_shp$POPMEAN_1

mwi_shp<-readOGR(dsn = path.expand("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/data/Malawi"),
                   layer = "mal_dist_join", stringsAsFactors = FALSE)
summary(mwi_shp)

rwa_shp<-readOGR(dsn = path.expand("C:/Users/Hari/Documents/PIH/consulting/Rwanda"),
                   layer = "RWA_export_merge", stringsAsFactors = FALSE)
summary(rwa_shp)

tz_shp<-readOGR(dsn = path.expand("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/TZ paper final/data"),
                layer = "tz_dist_join_v9", stringsAsFactors = FALSE)
summary(tz_shp)


tzdat<-tz_shp@data
summary(tzdat$T2MEAN)
tzdat<-tzdat[tzdat$T2MEAN>0,]


cor.test(tzdat$POPMEAN, tzdat$T2MEAN, method="spearman")
cor.test(tzdat$LOGPOPMEAN, tzdat$T2MEAN)

tz_shp@data$clust_HH_2 <-tz_shp@data$t2_fdr_HH
tz_shp@data$clust_LL_2 <-tz_shp@data$t2_fdr_LL
tz_shp@data$clust_LH_2 <-tz_shp@data$t2_fdr_LH
tz_shp@data$clust_HL_2 <-tz_shp@data$t2_fdr_HL



########################################
## Step 2: some functions
########################################

## need to create factor variable using cluster

getCluster<-function(indat){
  indat$cluster<-0
  indat$cluster[indat$clust_HH_2==1]<-1
  indat$cluster[indat$clust_LL_2==1]<-2
  indat$cluster[indat$clust_LH_2==1]<-3
  indat$cluster[indat$clust_HL_2==1]<-4
  
  indat$cluster.f<-factor(indat$cluster, levels = c("0", "1", "2","3","4"))
  ## 0 = NS, 1 = HH, 2 = LL, 3=LH, 4=HL
  return(indat)
    
}




plot_prep<-function(indat, cname){
  indat$country<-cname
  indat<-indat[,c("POPMEAN","T2MEAN","cluster.f","country")]
  return(indat)
}

## Process data files

ken_shp@data<-getCluster(ken_shp@data)
ken_plot<-plot_prep(ken_shp@data,"KEN")

mwi_shp@data<-getCluster(mwi_shp@data)
mwi_plot<-plot_prep(mwi_shp@data,"MWI")

rwa_shp@data<-getCluster(rwa_shp@data)
rwa_plot<-plot_prep(rwa_shp@data,"RWA")

tz_shp@data<-getCluster(tz_shp@data)
tz_plot<-plot_prep(tz_shp@data,"TZA")

## plots

go_small<-rbind(ken_plot,tz_plot,rwa_plot,mwi_plot)


scaleFUN <- function(x) sprintf("%.2f", x)


go_small$country <- factor(go_small$country, levels = c("KEN", "TZA", "RWA", "MWI"))

## need to create the separate population density and travel time means for each country
library(tidyverse)
go_ax <- go_small %>% 
  group_by(country) %>%
  summarise(tt = median(T2MEAN, na.rm=TRUE))

go_pop <- go_small %>%
  group_by(country) %>%
  summarise(pd = median(POPMEAN, na.rm=TRUE))


## Plots accounting for spatial structure
## change color scheme
## updated colors gray78 (NS) #440154FF (High/Long) #55C667FF  (Low/Long) 
#39568CFF (Low/Short) #FDE725FF (High/Short)
colvec<-c("gray78","#440154FF","#55C667FF","#39568CFF","#FDE725FF")

p1<-ggplot(go_small, aes(x = POPMEAN, y = T2MEAN, color=cluster.f)) + 
  geom_point() + 
  scale_x_continuous(trans = 'log', labels=scaleFUN) +
  scale_y_continuous(trans = 'log', labels=scaleFUN) +
  scale_color_manual(labels=c("NS", "High/Long","Low/Long","High/Short"), values = c("gray78","#440154FF", "#39568CFF","#FDE725FF")) +
  # colors are messed up, I have done this manually
  labs(x=expression('Average District Population per 10,000m'^2), y="Time to Health Center (min)", labs(fill='Cluster')) +
  guides(color=guide_legend(title="Legend")) +
  geom_hline(data=go_ax, aes(yintercept=tt), linetype="dashed") +
  geom_vline(data=go_pop, aes(xintercept=pd), linetype="dashed") +
  geom_hline(yintercept = 120, colour="blue", linetype="solid") +
  theme_classic() +
  theme(legend.position = "bottom") + 
  facet_wrap(~country, nrow=2)

