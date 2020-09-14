#########################################################
# Concept Paper
# Get TZA, MWI, RWA, KEN
# Keep only 
# - sites that seem to have research/med
# - geocodes
#########################################################

library(tidyverse)

godat<-read.csv("C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/Concept Paper/data/go_east africa_geocode_clean.csv", stringsAsFactors = FALSE, header=TRUE)

getCountry<-function(indat,cstring){
  indat<-indat[indat$country_names==cstring,]
  indat<-indat[!is.na(indat$Latitude),]
  return(indat)
}

write.csv(getCountry(godat,"Kenya"),"C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/Concept Paper/data/research programs/ken_progs.csv", row.names = FALSE) # drop 1
write.csv(getCountry(godat,"Malawi"),"C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/Concept Paper/data/research programs/mwi_progs.csv", row.names = FALSE) # drop 0
write.csv(getCountry(godat,"Rwanda"),"C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/Concept Paper/data/research programs/rwa_progs.csv", row.names = FALSE) # drop 1
write.csv(getCountry(godat,"Tanzania"),"C:/Users/Hari/Documents/Harvard/Global Cancer Epidemiology/John Flanigan - Diagnostics/Concept Paper/data/research programs/tza_progs.csv", row.names = FALSE)


