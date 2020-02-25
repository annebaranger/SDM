####### this code is made to extract the occurence data from BIEN, ALA and GBIF ##########

#### specify a pathway in which you want to save the data ######
saving_folder="C:/Users/HP/Desktop/"



###### uploading data obtained from script 2
Final_Cosmetic_list=read.csv2("Final_cosmetic_list.csv",stringsAsFactors = FALSE)

setwd(saving_folder)

library(BIEN)
library(rgbif)
library(ALA4R)
library(plyr)

##### we create data frames with only one column, containing BIEN, GBIF or ALA  ######
######     in order to add a column: EXtracted_from, in the occurence data frame   #######
Bien<-data.frame(Extracted_from=as.character("BIEN"))
Gbif<-data.frame(Extracted_from=as.character("GBIF"))
Ala<-data.frame(Extracted_from=as.character("ALA"))




#####   with this loop, we extracted and save all the occurence data species per species  ####

for (i in 1:length(Final_Cosmetic_list) ){
  
  species=Final_Cosmetic_list$Name[i]
  restriction<-data.frame(restriction=Final_Cosmetic_list$restriction[i]) ### keep in mind the restriction applied to our species as ingredients
  Species_TNRS<-data.frame(Species_TNRS=species) ### putting species name in the a data frame that will represent the first column of the final table
  
  ### extracting occcurence data from BIEN
  occ_bien1=as.data.frame(BIEN_occurrence_species(species, cultivated = TRUE, only.new.world = FALSE, all.taxonomy = TRUE, native.status = TRUE, natives.only = FALSE, observation.type = TRUE, political.boundaries = TRUE, collection.info = TRUE))
  occ_bien=cbind(occ_bien1,Bien,Species_TNRS,restriction)
  names(occ_bien)[names(occ_bien) == 'scrubbed_species_binomial'] <- 'Species'
  names(occ_bien)[names(occ_bien) == 'locality'] <- 'Locality'
  names(occ_bien)[names(occ_bien) == 'country'] <- 'Country'
  names(occ_bien)[names(occ_bien) == 'latitude'] <- 'Latitude'
  names(occ_bien)[names(occ_bien) == 'longitude'] <- 'Longitude'
  names(occ_bien)[names(occ_bien) == 'date_collected'] <- 'Event_date'

  ### extracting occurence date from GBIF
  occ_gbif1=occ_search(scientificName = species)
  occ_gbif_data=as.data.frame(occ_gbif1$data)
  occ_gbif=cbind(occ_gbif_data,Gbif,Species_TNRS,restriction)
  names(occ_gbif)[names(occ_gbif) == 'species'] <- 'Species'
  names(occ_gbif)[names(occ_gbif) == 'locality'] <- 'Locality'
  names(occ_gbif)[names(occ_gbif) == 'country'] <- 'Country'
  names(occ_gbif)[names(occ_gbif) == 'decimalLatitude'] <- 'Latitude'
  names(occ_gbif)[names(occ_gbif) == 'decimalLongitude'] <- 'Longitude'
  names(occ_gbif)[names(occ_gbif) == 'eventDate'] <- 'Event_date'
  names(occ_gbif)[names(occ_gbif) == 'elevation'] <- 'Elevation'

  ### extracting occurence date from ALA
  occ_ala1=occurrences(species, download_reason_id =11)
  occ_ala_data=occ_ala1$data
  occ_ala=cbind(occ_ala_data,Ala,Species_TNRS,restriction)
  names(occ_ala)[names(occ_ala) == 'species'] <- 'Species'
  names(occ_ala)[names(occ_ala) == 'locality'] <- 'Locality'
  names(occ_ala)[names(occ_ala) == 'country'] <- 'Country'
  names(occ_ala)[names(occ_ala) == 'latitudeOriginal'] <- 'Latitude'
  names(occ_ala)[names(occ_ala) == 'longitudeOriginal'] <- 'Longitude'
  names(occ_ala)[names(occ_ala) == 'eventDate'] <- 'Event_date'
  names(occ_ala)[names(occ_ala) == 'maximumElevationInMeters'] <- 'Elevation'
  
  ### putting all data in a single table
  occ_tot=rbind.fill(occ_ala,occ_gbif,occ_bien)
  occ_tot_reduced=as.data.frame(cbind(occ_tot$Species_TNRS,
                                      occ_tot$Species,
                                      occ_tot$Longitude,
                                      occ_tot$Latitude,
                                      occ_tot$Country,
                                      occ_tot$Locality,
                                      occ_tot$Elevation,
                                      occ_tot$Event_date,
                                      as.data.frame(occ_tot$Extracted_from),
                                      occ_tot$restriction))
  
  names(occ_tot_reduced) <- c("Species_TNRS",
                              "Species",
                              "Longitude",
                              "Latitude",
                              "Country",
                              "Locality",
                              "Elevation",
                              "Event_date",
                              "Extracted_from",
                              "Restriction")
  
  #### occ_tot contains all the information we have on the species and the occurence records
  #### occ_tot_reduced only contains useful informations for the next steps
  write.csv2(occ_tot, paste("occ_tot_",as.character(Final_Cosmetic_list[i,1]),".csv"), row.names=FALSE, sep="t",dec=",", na=" ")
  write.csv2(occ_tot_reduced, paste("occ_reduced_",as.character(Final_Cosmetic_list[i,1]),".csv"), row.names=FALSE, sep="t",dec=",", na=" ")
  }
  