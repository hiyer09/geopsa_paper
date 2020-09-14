# Code and Data for Geographic-Population Services Access Model Paper

This folder contains datasets and R code files to produce the scatter plots and maps in our paper. Below are some details regarding how the files are organized.

1. geo_psa/raster 

This folder contains travel time rasters produced using Access Mod 5, organized by country using their three-letter abbreviations. Within each folder, you will find rasters corresponding to travel time to closest health center and cancer referral center. For more details regarding the methods, please refer to the paper.

2. geo_psa/code

This folder contains R code used to generate the bivariate local indicator of spatial autocorrelation scatter plots. These code files read in spatial datasets that were processed in ArcMap (zonal statistics) and GeoDa (bivariate LISA) and can be adapted to data files in your specific context.
