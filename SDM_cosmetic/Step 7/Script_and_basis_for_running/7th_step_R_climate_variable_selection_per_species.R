library(sf)
library(raster)
library(rgdal)
library(tidyverse)
library(rgeos)
library(scales)
library(fasterize)
library(maptools)
library(biomod2)
library(usdm)
library(plyr)

saving_folder="C:/Users/190384/Desktop"

### where is your species list
species_list_path="C:/Users/190384/Desktop"
setwd(species_list_path)

Species_list=read.csv2("Final_cosmetic_list.csv",stringsAsFactors = FALSE)
Species_list=as.data.frame(Species_list[order(Species_list[,1],decreasing=F),])


###### keeping the climate variables that are not too much correlated

### where are the raw climate data from CHelsa_1979_2013
path="S:/BD_SIG/climat/monde/Chelsa_1979_2013/climatologies/bio"

files <- list.files(path,pattern="tif", full.names=TRUE )

###getting rid of unused files/climate variables
files <- files[-c(12:13)]
predictors <- stack(files)

### where are your corrected occurence data
occ_corrected_folder="C:/Users/190384/Desktop/Occurrences_corrected"
setwd(occ_corrected_folder)

### creation of an empty data frame that will be filled with the selected climate variables for each species
Predictor_list=structure(list(character()), class = "data.frame")

for (k in 1:nrow(Species_list)) {
  species=Species_list[k,1]
  setwd(occ_corrected_folder)
  Species_occ=read.csv2(paste("occ_corrected_", species,".csv"),stringsAsFactors = FALSE)
  SpXY=Species_occ[,1:2]
  
  ### if there is at least 10 occurrence records for species i ...
  if ((nrow(SpXY)>9)==TRUE){
    Sp_predictors<-raster::extract(predictors,SpXY) ### extraction of climate variables values at occurence points
    Sp_predictors=as.data.frame(Sp_predictors)
    Sp_predictors_filtered=as.data.frame(na.omit(Sp_predictors)) ### getting rid of potentila NA
    Sp_predictors_cor_test=vifcor(Sp_predictors_filtered,th=0.7) ### correlation test between climate variables and selection
                                                                 ### of the right subset of variables (not too much correlated) 
    Sp_predictor_used=as.data.frame(Sp_predictors_cor_test@results) ### results of the test
    Sp_predictor_used_list=as.data.frame(Sp_predictor_used[,-2]) ### only keeping variable names
    colnames(Sp_predictor_used_list)="Predictor_variables"
    Sp_predictors_list=as.data.frame(cbind(species,t(Sp_predictor_used_list))) ### adding species name in the first column, before variable names
    
    Predictor_list=as.data.frame(rbind.fill(Predictor_list,Sp_predictors_list))} ### For each species, 
                                                                                 ### list of variables to use as predictors
  
  print(k)
}


setwd(saving_folder)
write.csv2(Predictor_list,"predictor_list.csv",row.names = F,sep= ";")
