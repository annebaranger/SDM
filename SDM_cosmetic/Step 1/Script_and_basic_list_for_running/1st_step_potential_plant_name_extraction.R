##### First code to put the initial cosing data in the right format in order to send it to tnrs ####
#####   for comparison   ######

### name of a folder in which your final lists will be saved 
saving_folder="C:/Users/HP/Desktop"

####### uploading the initial dataset : list of 26109 cosmetic ingredients ########
####### (cosing-ingredients-frangrance inventory, last update 01/01/2019)  ########

cosmetic_list1=read.csv2("Cosmetic_list.csv")
cosmetic_list1=as.data.frame(cosmetic_list1)
cosmetic_list1 <- cosmetic_list1[ , - c(3:6,10:11)]       ### getting rid of useless columns
colnames(cosmetic_list1) = c('COSING Ref No', "INCI Name", "Chem/IUPAC Name/Description", "Restriction", "Function")     ### getting rid of useless rows


####### data treatment to take the potential species names out of the list of ingredients  ###########
#######                 through pattern recognition                       #########

library(stringr)
L<- data.frame(number= integer(0),name= numeric(0),description=character(0), restriction=character(0), fonction=character(0))  #### creation of an empty data frame to fill with the extracted words (column: nom) and the description that goes with #######
for (i in 5:nrow(cosmetic_list1)){            ##### for all rows in the data frame
  x=cosmetic_list1[i,2]        ##### we will study the second column where we can find species names
  a=str_detect(x, "()")        ##### a contains "TRUE" if there are parentheses in the box (may signify a species name is in the box)
  b=str_detect(x, "/")         ##### b contains "TRUE" if there are slashes in the box (may signify a species name is in the box)
  if ((a==FALSE) | (b==FALSE)){       ###### if a and b are "FALSE" (no parenthesis nor slash) the species name id in the two first words (if there is one)
    x1=as.data.frame(str_split(as.character(x),  ' ' ))   ##### we split the 'sentence' in the box in a data frame in which each boxes is a word
    if (length(x1[,1])>1) {       ###### if there are more than two words in the 'sentence', we only keep the first two (probable genus and species name)
      x2=as.data.frame(cbind(as.integer(cosmetic_list1[i,1]),as.character(paste(x1[1,1],x1[2,1])),as.character(cosmetic_list1[i,3]),as.character(cosmetic_list1[i,4]),as.character(cosmetic_list1[i,5])))  ##### we create a data frane (x2) that contains all the information on the potential species (restriction, function, ...) and the potential species name (=the two first words (= the two first raws of the x1 data frame)) in the second columns
      colnames(x2)<-cbind("number", "name", "description", "restriction", "fonction")    ##### we rename the colums
      L=rbind(L,x2)}    ###### we add the raw that contains the potential species name to the new data frame
    else {     ####### if there are less than two words in the 'sentence', we keep the only words as a potential species or genus name
      x2=as.data.frame(cbind(as.integer(cosmetic_list1[i,1]),x1,as.character(cosmetic_list1[i,3]),as.character(cosmetic_list1[i,4]),as.character(cosmetic_list1[i,5])))    ##### we create a data frane (x2) that contains all the information on the potential species (restriction, function, ...) and the potential species name in the second columns
      colnames(x2)<-cbind("number", "name", "description", "restriction", "fonction")  ##### we rename the columns
      L=rbind(L,x2)}}   ###### we add the raw that contains the potential species name to the new data frame
  else {       #######  if a OR b are "TRUE" (there are slashes or parentheses)
    x1=as.data.frame(str_split(as.character(x),  '/' ))    ##### we firstly separate the sentence with the slashes (that separates species names if there are several ones). And we create a data frame,x1, with the various 'series of words' it extracted
    for (j in 1:length(x1[,1])){   ####### for all the 'series of words' containing a potential species name
      x2=as.data.frame(str_split(as.character(x1[j,1]),  ' ' ))    ##### we separates the words and create a data frame (x2) with one word per box
      if (length(x2[,1])>1) {    ###### if there are more than two words in the 'series of words',    
        x3=as.data.frame(cbind(as.integer(cosmetic_list1[i,1]),as.character(paste(x2[1,1],x2[2,1])),as.character(cosmetic_list1[i,3]),as.character(cosmetic_list1[i,4]),as.character(cosmetic_list1[i,5])))    ##### we only keep the first two as potenetial species name and we create a raws that will be added to the new data frame with all the information needed
        colnames(x3)<-cbind("number", "name", "description", "restriction", "fonction")   ##### we rename the columns
        L=rbind(L,x3)}     ###### we add the raw that contains the potential species name to the new data frame
      else {     ###### if there are less than two words, we take the only words as a potential species or genus name
        x3=as.data.frame(cbind(as.integer(cosmetic_list1[i,1]),x2,as.character(cosmetic_list1[i,3]),as.character(cosmetic_list1[i,4]),as.character(cosmetic_list1[i,5])))
        colnames(x3)<-cbind("number", "name", "description", "restriction", "fonction")    ##### we rename the columns
        L=rbind(L,x3)    ###### we add the raw that contains the potential species name to the new data frame
      } }}}

####### At the end, L is a list of 33381 series of two words containing, for some, the name of a plant species #############


##### To save the list on the desktop ######

cometic_2_words_list=L   ###### cometic_2_words_list is the L list (33381 series of two words) #####
setwd(saving_folder)
write.csv2(cometic_2_words_list, "2_words_list.csv", row.names=FALSE, sep="t",dec=",", na=" ")

##### We get rid of the parentheses that could be attached to some of the words and could be an issue for tnrs comparison with the plant species list #####
##### This is donne through Excel and the new list obtained is called cosmetic_2_words_list_without_parenthesis and is uploaded from a csv file #####
cosmetic_2_words_list_without_parenthesis=read.csv2("2_words_list_without_parenthesis.csv")
library(dplyr)

##### we do a unique to be sure we have no repetition of the same serie of two words #####
cosmetic_2_words_list_wp_unique=as.data.frame(cosmetic_2_words_list_without_parenthesis %>% distinct(cosmetic_2_words_list_without_parenthesis$name, .keep_all = TRUE))  ### J is a data frame with no repetition of potential species, there are 19102 different series of two words
write.csv2(cosmetic_2_words_list_wp_unique, "2_words_list_wp__unique.csv", row.names=FALSE, sep="t",dec=",", na=" ")


#######  we organise the data to send them to tnrs for comparison, we only keep the columns with the series of two words ####
list_for_tnrs=cosmetic_2_words_list_wp_unique$name
write.csv2(list_for_tnrs, "list_for_tnrs.csv", row.names=FALSE, sep="t",dec=",", na=" ")











