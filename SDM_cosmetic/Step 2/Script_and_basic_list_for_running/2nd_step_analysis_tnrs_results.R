#######  Second code to analyse the data received from tnrs   #########$

saving_folder="C:/Users/HP/Desktop"


##### Uploading the data from tnrs : only best matches ######
tnrs_best1=read.csv2("tnrs_results_best_matches.csv",stringsAsFactors = FALSE)   
tnrs_best1=tnrs_best1[-1,]

##### Uploading the data from obtained from script 1 ######
cosmetic_2_words_list_wp_unique=read.csv2("2_words_list_wp__unique.csv",stringsAsFactors = FALSE)

##### adding the information on restriction, function and description for all species #####
tnrs_best1_all_informations=cbind(tnrs_best1,cosmetic_2_words_list_wp_unique[,3:5])
setwd(saving_folder)
write.csv2(tnrs_best1_all_informations, "tnrs_results_information_on_cosmetics.csv", row.names=FALSE, sep="t",dec=",", na=" ")

##### In order to get rid of useless information or ingredients that are not species, we choose to 
##### discriminate the results with an overall score superior to 0.95 and,
##### for all the scores between 0.95 and 0.80, we checked for potential species by hand, directly in 
##### the excel file ######
tnrs_best_sup_0.95=subset(tnrs_best1_all_informations,Overall_score>0.95)   ### we only keeep the best matches with an overall score superior to 0.95
tnrs_best_inf_0.95_ok=read.csv2("list_inf_0.95_but_still_good.csv", stringsAsFactors=FALSE)   ### we upload the data hand cheked on the excel file
TNRS_global_list_species=rbind(tnrs_best_sup_0.95,tnrs_best_inf_0.95_ok)  ##### we concatenate both list to have the global list of species used in cosmetics


##### Organising the data to have a well designed table ####

### Keeping the accepted name for each species ###
A<-data.frame(Name=as.character(0)) #### we create a data frame to add the accepted name, it will be a new column for the last TNRS_global_list_species 
TNRS_global_list_species_2=cbind(TNRS_global_list_species,A)  ### we add the new empty column to the former data frame

##### We will do the difference between taxonomic status #####
####   If it is "accepted", "synonym", or "misapplied", we use the accepted name as new global name : 2894 raws ####

TNRS_species_with_accepted_names=subset(TNRS_global_list_species_2,Taxonomic_status=="Accepted"|Taxonomic_status=="Synonym"|Taxonomic_status=="misapplied")   #### we only keep the species for which there is an accepted name or a synonym or misapplication ####
TNRS_species_with_accepted_names$Name=TNRS_species_with_accepted_names$Accepted_name #### and we put the accepted name in the column name ####


#### If it is "rejected", there is "no opinion", it is "invalid", "illegitimate" or else ...227 raws #####

TNRS_matched_names=subset(TNRS_global_list_species_2,Taxonomic_status=="Invalid"|Taxonomic_status==""|Taxonomic_status=="Illegitimate"|Taxonomic_status=="No opinion"|Taxonomic_status=="Rejected name")  ### we only keep the species for which there is a taxonomic status corresponding to the conditions ####
#### Among these, we do a subset for those wo have an Accepted name : 15 raws and we keep this accepted name as the new one ### ####
TNRS_matched_accepted_names=subset(TNRS_matched_names,Accepted_name!="")  
TNRS_matched_accepted_names$Name=TNRS_matched_accepted_names$Accepted_name
#### And for those wo have no accepted name, we use the name matched : 212 raws #####
TNRS_matched_matched_names=subset(TNRS_matched_names,Accepted_name=="")
TNRS_matched_matched_names$Name=TNRS_matched_matched_names$Name_matched

####  We combine all the previous lists to get the global list of the 3121 species : Cosmetic_list_with_names  ####
Cosmetic_list_with_names=rbind(TNRS_species_with_accepted_names,TNRS_matched_accepted_names,TNRS_matched_matched_names)
#### we only keep the interesting columns and put them in the right order  ####
Final_Cosmetic_list_with_names=Cosmetic_list_with_names[,c(40,29,32,14,16,31,26,3,34,39,38,37)]
#### we do a unique so that no species is repeated : 2887 raws dont 2579 species  #####
library(dplyr)
Final_Cosmetic_list_with_unique_names=as.data.frame(Final_Cosmetic_list_with_names %>% distinct(Name, .keep_all = TRUE))
write.csv2(Final_Cosmetic_list_with_unique_names, "Final_cosmetic_list_with_genus.csv", row.names=FALSE, sep="t",dec=",", na=" ")


####  we will use a list of all the species and get rid of the raws that contains only a genus ###
Final_Cosmetic_list=subset(Final_Cosmetic_list_with_unique_names, Accepted_name_rank=="species")
write.csv2(Final_Cosmetic_list, "Final_cosmetic_list_update.csv", row.names=FALSE, sep="t",dec=",", na=" ")



list_before=read.csv2("Final_cosmetic_list.csv")
list_before_species_unique=subset(list_before, Accepted_name_rank=="species")
write.csv2(list_before_species_unique, "Final_cosmetic_list_non_update.csv", row.names=FALSE, sep="t",dec=",", na=" ")

list1=read.csv2("Final_cosmetic_list_non_update.csv")
list2=read.csv2("Final_cosmetic_list_update.csv")

L=list2[!((list2$Name %in% list1$Name)), ]
write.csv2(L, "species_added_in_the_updated_version.csv", row.names=FALSE, sep="t",dec=",", na=" ")


