
####### R code to clean the occurrences dat ######
###### see Alexzander Zizka coordinateCleaner package  #####


saving_folder_summary_flagged_records="C:/Users/HP/Desktop"
saving_folder_occ_corrected="C:/Users/HP/Desktop"

### where are your reduced occurrence data obtained step 3 ####
occ_reduced_folder="C:/Users/HP/Desktop"

### where is your final cosmetic list ####
basic_list_folder="C:/Users/HP/Desktop/Liste de base/"

library(readxl)
library(CoordinateCleaner)


###### downloading the list with all the species names #######

setwd(basic_list_folder)
Final_Cosmetic_list=read.csv2("Final_cosmetic_list.csv",stringsAsFactors = FALSE)


###### Creating a data frame in which you will have the percentage of flagged records according to the test ######
######     for each species and the number of good records   ########
SUMMARY=data.frame("Species"=as.numeric(0),
                   "coordinate_validity"=as.numeric(0),
                   "conversion_errors"=as.numeric(0),
                   "zero_coordinates"=as.numeric(0),
                   "country_capitals"=as.numeric(0),
                   "country_centroids"=as.numeric(0),
                   "duplicates"=as.numeric(0),
                   "equal_lat_lon"=as.numeric(0),
                   "GBIF_headquarters_and_flagged_records_around_Copenhagen"=as.numeric(0),
                   "biodiversity_institutions"=as.numeric(0),
                   "sea_coordinates"=as.numeric(0),
                   "geographic_outliers"=as.numeric(0),
                   "percentage_of_flagged_records"=as.numeric(0),
                   "number_of_records_ok"=as.numeric(0))

###### Creating an empty data frame that will be filled in the for loop, each time you change species you add ###
######   a row to the data frame  #####
Final_summary<-NULL


#####   Loop to clean the occurrence dataset for each species, only keeping the good ones and the percentage  ####
##### of records flagged by each test   ########"
for (i in 1:length(Final_Cosmetic_list)){
  
  species=Final_Cosmetic_list$Name[i]  ##### name of the species i
  
  SUMMARY$Species=as.character(species) #### filling the data.frame created earlier to save the percentage of records flagged per test and species
  
  TEST_species=paste("occ_reduced_",species,".xlsx")   #### name of the file in which there are the occurrences 
                                                       #### of the species that will be tested
  
  #### name of the path in wich you saved your occurence data
  setwd(occ_reduced_folder)
  occ_TEST=read_excel(TEST_species)  ##### data_frame with the occurrences data to be tested
  
  ######  Removing non-numeric and not available coordinates  #######
  ######  as well as lat >90, la <-90, lon > 180 and lon < -180   #######
  
  TEST=cc_val(occ_TEST, lon = "Longitude", lat = "Latitude",
              value = "clean", verbose = TRUE)
  SUMMARY$coordinate_validity=(nrow(occ_TEST)-nrow(TEST))/nrow(occ_TEST)   #### give the percentage of records 
                                                                           #### flagged by this test
  
  
  ######  Removing erroneous conversion from a degree minute format to a decimal degree format  ######
  
  TEST1=clean_dataset(TEST, lon = "Longitude", lat = "Latitude",
                      ds = "Species_TNRS",tests = c("ddmm"),
                      value = "clean",verbose= TRUE)
  SUMMARY$conversion_errors=(nrow(TEST)-nrow(TEST1))/nrow(occ_TEST)
  
  
  ###### here, if there is 0 occurrences remaining we stop the testing and we set the following values #####
  
  if ((nrow(TEST1)==0)==TRUE){ 
    SUMMARY$zero_coordinates=0
    SUMMARY$conversion_errors=(nrow(TEST)-nrow(TEST1))/nrow(occ_TEST)  
    SUMMARY$country_capitals=0
    SUMMARY$country_centroids=0
    SUMMARY$duplicates=0
    SUMMARY$equal_lat_lon=0
    SUMMARY$GBIF_headquarters_and_flagged_records_around_Copenhagen=0
    SUMMARY$biodiversity_institutions=0
    SUMMARY$sea_coordinates=0
    SUMMARY$geographic_outliers=0
    SUMMARY$percentage_of_flagged_records=1
    SUMMARY$number_of_records_ok=0    ### there is no remaining records
    
    setwd(saving_folder_occ_corrected)
    write.csv2(TEST1, paste("occ_corrected_",as.character(species),".csv"), row.names=FALSE)
    
    }
  
  
  
  ##### but if there are still some occurrences in the dataset we continue the test  ######
  
  if ((nrow(TEST1)==0)==FALSE){  
    
  #####   to avoid problem with missing coordinates and NA, see Zizka mail ####  
  TEST1$Longitude=as.numeric(TEST1$Longitude)
  TEST1$Latitude=as.numeric(TEST1$Latitude)
  TEST2=TEST1[!is.na(TEST1$Longitude),]
  TEST3=TEST2[!is.na(TEST2$Latitude),]
  
  #####  Removing records with either zero longitude or latitude 
  #####   and a radius of 0.1 decimal degree around the point at zero longitude and zero latitude
  TEST4=cc_zero(TEST3, lon = "Longitude", lat = "Latitude",
                value = "clean", verbose = TRUE)
  SUMMARY$zero_coordinates=(nrow(TEST1)-nrow(TEST4))/nrow(occ_TEST)  #### give the percentage of records 
                                                                     #### flagged by this test
  
  
  ##### Removing records within a 10km radius around country capitals  ####
  TEST5=cc_cap(TEST4, lon = "Longitude", lat = "Latitude",
               species = "Species_TNRS", value = "clean", verbose = TRUE)
  SUMMARY$country_capitals=(nrow(TEST4)-nrow(TEST5))/nrow(occ_TEST)   #### give the percentage of records 
                                                                      #### flagged by this test
  
  
  ##### Removing records within a 1 kmradius around the geographic centroids of political countries and provinces ###
  TEST6=cc_cen(TEST5, lon = "Longitude", lat = "Latitude",
               species = "Species_TNRS",test = "both", value = "clean", verbose = TRUE)
  SUMMARY$country_centroids=(nrow(TEST5)-nrow(TEST6))/nrow(occ_TEST)  #### give the percentage of records 
                                                                      #### flagged by this test
  
  ##### Removing duplicated records #####
  TEST7=cc_dupl(TEST6, lon = "Longitude", lat = "Latitude",
                species = "Species_TNRS", value = "clean", verbose = TRUE)
  SUMMARY$duplicates=(nrow(TEST6)-nrow(TEST7))/nrow(occ_TEST)      #### give the percentage of records 
                                                                   #### flagged by this test
  
  ##### Removing records with equal latitude and longitude coordinates, either exact or absolute ####
  TEST8=cc_equ(TEST7, lon = "Longitude", lat = "Latitude",
               test = "absolute", value = "clean", verbose = TRUE)
  SUMMARY$equal_lat_lon=(nrow(TEST7)-nrow(TEST8))/nrow(occ_TEST)  #### give the percentage of records 
                                                                  #### flagged by this test
  
  ##### records within 0.5 degree radius around the GBIF headquarters in Copenhagen  #####
  TEST9=cc_gbif(TEST8, lon = "Longitude", lat = "Latitude",
                species = "Species_TNRS", value = "clean", verbose = TRUE)
  SUMMARY$GBIF_headquarters_and_flagged_records_around_Copenhagen=(nrow(TEST8)-nrow(TEST9))/nrow(occ_TEST)  #### give the percentage of records 
                                                                                                            #### flagged by this test
  
  ##### Removing records assigned to the location of zoos, botanical gardens, herbaria, universities and museums ###
  TEST10=cc_inst(TEST9, lon = "Longitude", lat = "Latitude",
                 species = "Species_TNRS", value = "clean", verbose = TRUE)
  SUMMARY$biodiversity_institutions=(nrow(TEST9)-nrow(TEST10))/nrow(occ_TEST)   #### give the percentage of records 
                                                                                #### flagged by this test
  
  ##### Removing coordinates outside the reference landmass #####
  TEST11=cc_sea(TEST10, lon = "Longitude", lat = "Latitude",
                ref = NULL, value = "clean", verbose = TRUE)
  SUMMARY$sea_coordinates=(nrow(TEST10)-nrow(TEST11))/nrow(occ_TEST)   #### give the percentage of records 
                                                                       #### flagged by this test
  
  ##### Removing records that are outliers in geographic space (erroneous coordinates, ######
  ######  for example due to data entry errors, imprecise geo-references, individuals in horticulture/captivity) ####
  TEST12=cc_outl(TEST11, lon = "Longitude", lat = "Latitude",
                 species = "Species_TNRS",value="clean")
  SUMMARY$geographic_outliers=(nrow(TEST11)-nrow(TEST12))/nrow(occ_TEST) #### give the percentage of records 
                                                                         #### flagged by this test
  
  #### give the total percentage of records flagged ######
  SUMMARY$percentage_of_flagged_records=(nrow(occ_TEST)-nrow(TEST12))/nrow(occ_TEST) 
  
  #### give the number of good occurrences ######
  SUMMARY$number_of_records_ok=nrow(TEST12)
  
  #####  saving the file with the occurrences corrected  #####
  setwd(saving_folder_occ_corrected)
  write.csv2(TEST12, paste("occ_corrected_",as.character(species),".csv"), row.names=FALSE)
  
  }
  
  ####  adding the results of species i to the summary of fagged records #####
  Final_summary=rbind(Final_summary,SUMMARY)
  print(i)

}

######  saving the summary on the flagged records ######
setwd(saving_folder_summary_flagged_records)
write.csv2(Final_summary, "Final_Summary.csv", row.names=FALSE)
