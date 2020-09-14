#################################################################
# R Code for Bivariate LiSA Scatter - Concept Paper
# Objective: Evaluate urban-rural access to cancer research centers
# in SSA
# Goals: read in all GeoDA shapefiles
# Do the combined scatter plot, with factor variable for different
# cluster types


# Updated 6/11/20 - restricted to cancer referral sites, 
# applied FDR (B-H correction)
# Change the axis to reflect in-country - updated

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

tz_shp@data$GOMEAN2<-tz_shp@data$GOMEAN


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
  indat<-indat[,c("POPMEAN","GOMEAN2","cluster.f","country")]
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
  summarise(tt = median(GOMEAN2, na.rm=TRUE))

go_pop <- go_small %>%
  group_by(country) %>%
  summarise(pd = median(POPMEAN, na.rm=TRUE))


## Plots accounting for spatial structure
## updated colors gray78 (NS) #440154FF (High/Long) #55C667FF  (Low/Long) 
#39568CFF (Low/Short) #FDE725FF (High/Short)
colvec<-c("gray78","#440154FF","#55C667FF","#39568CFF","#FDE725FF")



p3<-ggplot(go_small, aes(x = POPMEAN, y = GOMEAN2, color=cluster.f)) + 
  geom_point() + 
  scale_x_continuous(trans = 'log', labels=scaleFUN) +
  scale_y_continuous(trans = 'log', labels=scaleFUN) +
  scale_color_manual(labels=c("NS", "High/Long","Low/Short","Low/Long","High/Short"), values = colvec) +
  # colors are messed up, I have done this manually
  labs(x=expression('Average District Population per 10,000m'^2), y="Time to Cancer Referral Center (min)", labs(fill='Cluster')) +
  guides(color=guide_legend(title="Legend")) +
  geom_hline(data=go_ax, aes(yintercept=tt), linetype="dashed") +
  geom_vline(data=go_pop, aes(xintercept=pd), linetype="dashed") +
  geom_hline(yintercept = 120, colour="blue", linetype="solid") +
  theme_classic() +
  theme(legend.position = "bottom") + 
  facet_wrap(~country, nrow=2)

p3
## what's up with TZ

tzdat<-tz_shp@data
summary(tzdat$GOMEAN2)
tzdat<-tzdat[tzdat$GOMEAN2>0,]


cor.test(tzdat$POPMEAN, tzdat$GOMEAN2, method="spearman")
cor.test(tzdat$LOGPOPMEAN, tzdat$LOGGOMEAN2)
## there are lakes in the analysis.. they are probably already being dropped from
## GeoDa, but should confirm
