library(raster)

#### where are your climate data obtained from previous scripts
climate_data_path="C:/Users/190384/Desktop/climate_aggreg/"
### where is your species list
species_list_path="C:/Users/190384/Desktop/"
### where are your occurences data from step 4 
occurrences_corrected_path="C:/Users/190384/Desktop/"
saving_folder_occ_aggregated="C:/Users/190384/Desktop/"

setwd(climate_data_path)
### dowloading one raster containing climate data (just to have a  raster with the right cell size)
climat=raster(paste0(climate_data_path,"CHELSA_bio10_4_aggregate5.tif"))

setwd(species_list_path)
### downloading species_list 
species_list=read.csv2(paste0(species_list_path,"Final_cosmetic_list.csv"),stringsAsFactors = FALSE)


###, aggregatin the occurrence in order to have a maximum of one occurrence point per raster cell
for (i in 1:length(species_list)){
  species=species_list$Name[i]
  species_occ=read.csv2(paste0(occurrences_corrected_path,"occ_corrected_ ",species," .csv"),stringsAsFactors = FALSE)
  species_XY=species_occ[,1:2]
  species_XY_cell=cbind(as.data.frame(cellFromXY(climat,species_XY)),species_XY)
  species_XY_cell=species_XY_cell[! duplicated(species_XY_cell[, 1]), ]
  species_XY_aggregate=as.data.frame(species_XY_cell[,2:3])
  setwd(saving_folder_occ_aggregated)
  write.csv2(species_XY_aggregate,paste0("occ_aggregated_",species,".csv"),row.names=FALSE)
  print(i)
  }


