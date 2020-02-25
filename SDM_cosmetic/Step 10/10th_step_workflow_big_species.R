library(sf)
library(raster)
library(rgdal)
library(tidyverse)
library(rgeos)
library(scales)
library(fasterize)
library(maptools)
library(biomod2)
library(dplyr)
library(dismo)
library(ecospat)
library(stringr)
#library(rJava)

Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre1.8.0_241')
Sys.setenv(PATH='C:/Program Files/Java/jre1.8.0_241/bin/')

### parent foler in which the results of all your models will be saved
parent_folder='D:/Users/190384/Documents/Model/modele_report/'
setwd(parent_folder)


### where have you saved your climate rasters obtained step 7
raster_path="D:/Write/Noirault/Data/climate_extraction_per_spe/"

### where are your future climate rasters
future_climate_path="D:/Write/Noirault/Data/climat_tif/"


### where have you saved your predictors_list obtained step 9
predictor_list_path="D:/Users/190384/Documents/Model/"

### where is your list of variable filenames (not for the first time you run the script, see below)
predictor_filename_list_path="F:/Adrien/"

### where are your aggregated occurences data obtained step 8
occ_aggregated_path="D:/Write/Noirault/Data/occ_aggregated/"

### where is maxent.jar
path_to_maxent.jar='D:/Users/190384/Documents/Model/modele_report/'

### for species with errors
sp_with_errors_folder="D:/Users/190384/Documents/Model/modele_report/Species_with_errors/"


############### used fonction  #####################
####################################################
####################################################

##########################################################################################
### create a mask given the range of climate conditions in which you find the species ###
### the model will be applied to a region corresponding to this mask ###
### the rule is + or - 1 for standard deviation of each variable ###
######################################################################

extrapolationMask=function(thisEnv,
                           dat=NULL,
                           rule='sd',
                           val=1,
                           printReport=FALSE
){## for testing 					 
  # dat=bg.dat; thisEnv=env.new; rule='sd'; val=1
  #out=try({
  if(!rule=='sd') {
    stop("Only rule='sd' is implemented")
  } else {
    vars=colnames(dat)
    check=names(thisEnv) %in% vars
    if(any(!check)) warning(paste0('Some variables were not masked because they were not included in dat: ',names(thisEnv)[!check])) 
    masks=raster::stack(lapply(vars,function(x){
      e.min=min(dat[,x],na.rm=T)
      e.max=max(dat[,x],na.rm=T)
      e.sd=sd(dat[,x],na.rm=T)
      e.min.lim=e.min-val*e.sd
      e.max.lim=e.max+val*e.sd
      mask1=thisEnv[[x]] >= e.min.lim & thisEnv[[x]] <= e.max.lim
      raster::values(mask1)[raster::values(mask1)==0]=NA
      mask1
    }))
    names(masks)=vars
    #masks=raster::trim(masks)
    maskedEnv=raster::stack(lapply(vars,function(x){
      raster::mask(thisEnv[[x]],masks[[x]],maskValue=0)
    }))
    names(maskedEnv)=vars
  }
  
  # make a report of which layers lead to masking
  if(printReport){
    nonNAPreMask=apply(raster::values(thisEnv),2,function(x) sum(!is.na(x)))
    nonNAMasked=apply(raster::values(maskedEnv),2,function(x) sum(!is.na(x)))
    print(data.frame(nonNAPreMask=nonNAPreMask,nonNAMasked=nonNAMasked, numCellsMasked=nonNAPreMask-nonNAMasked))
  }
  return(list(masks=masks,maskedEnv=maskedEnv))
  #})
  #return(out)
}



###############################################################
##### Biomod_ensemble forecasting ############################
##############################################################


BIOMOD_EnsembleForecasting=function (EM.output=NULL,
                                     projection.output = NULL,
                                     new.env = NULL, 
                                     xy.new.env = NULL,
                                     proj.name = NULL, 
                                     binary.meth = NULL,
                                     filtered.meth = NULL,
                                     compress = NULL) 
{
  # .bmCat("Do Ensemble Models Projections")
  # args <- list(...)
  # args_checked <- .BIOMOD_EnsembleForecasting.check.args(EM.output, 
  #                                                        projection.output, new.env, selected.models, proj.name, 
  #                                                        total.consensus = FALSE, binary.meth, filtered.meth)
  # proj.name <- args_checked$proj.name
  # selected.models <- args_checked$selected.models
  
  ## add
  selected.models = get_built_models(EM.output)
  
  
  # output.format <- args$output.format
  
  ## add
  output.format = NULL
  
  # compress <- args$compress
  # do.stack <- args$do.stack
  # keep.in.memory <- args$keep.in.memory
  # on_0_1000 <- args$on_0_1000
  
  ## add
  do.stack=NULL
  keep.in.memory=NULL
  on_0_1000=NULL
  
  
  if (is.null(output.format)) {
    if (length(projection.output)) {
      if (projection.output@type != "RasterStack") 
        output.format <- ".RData"
      else output.format <- ".grd"
    }
    else {
      if (!inherits(new.env, "Raster")) 
        output.format <- ".RData"
      else output.format <- ".grd"
    }
  }
  if (is.null(compress)){ 
    compress <- FALSE}
  if (is.null(do.stack)) {
    do.stack <- TRUE
    if (!is.null(projection.output)) {
      if (all(grepl("individual_projections", projection.output@proj@link))) {
        do.stack <- FALSE
      }
    }
  }
  if (is.null(keep.in.memory)) {
    keep.in.memory <- TRUE
    if (!is.null(projection.output)) {
      keep.in.memory <- projection.output@proj@inMemory
    }
  }
  if (is.null(xy.new.env)) {
    if (!is.null(projection.output)) {
      xy.new.env <- projection.output@xy.coord
    }
    else {
      xy.new.env <- matrix()
    }
  }
  if (is.null(on_0_1000)){
    on_0_1000 = TRUE
    # rm(list = c("args_checked", "args"))
    proj_out <- new("BIOMOD.projection.out", proj.names = proj.name, 
                    sp.name = EM.output@sp.name, expl.var.names = EM.output@expl.var.names, 
                    models.projected = selected.models, xy.coord = xy.new.env, 
                    modeling.object.id = EM.output@modeling.id)
    proj_out@modeling.object@link = EM.output@link
    proj_is_raster <- FALSE}
  if (inherits(new.env, "Raster")) {
    proj_is_raster <- TRUE
  }
  # else if (length(projection.output)) {
  #   if (inherits(projection.output@proj, "BIOMOD.stored.raster.stack")) {
  #     proj_is_raster <- TRUE
  #   }
  # }
  if (proj_is_raster) {
    proj_out@proj <- new("BIOMOD.stored.raster.stack")
  }
  # else {
  #   proj_out@proj <- new("BIOMOD.stored.array")
  #   do.stack = TRUE
  # }
  dir.create(file.path(parent_folder,EM.output@sp.name, paste0("proj_", proj.name, 
                                                                          sep = "")), showWarnings = FALSE, recursive = TRUE, 
             mode = "777")
  indiv_proj_dir <- file.path(parent_folder,EM.output@sp.name, paste0("proj_", 
                                                                                 proj.name, sep = ""), "individual_projections")
  dir.create(indiv_proj_dir, showWarnings = FALSE, recursive = TRUE, 
             mode = "777")
  saved.files <- c()
  needed_predictions <- get_needed_models(EM.output, selected.models = selected.models)
  if (length(projection.output)) {
    formal_pred <- get_predictions(projection.output, full.name = needed_predictions,
                                   as.data.frame = ifelse(projection.output@type ==
                                                            "array", T, F))
  }
  else {
    tmp_dir <- paste0("Tmp", format(Sys.time(), "%s"), sep = "")
    formal_pred <- BIOMOD_Projection(modeling.output = load_stored_object(EM.output@models.out.obj), 
                                     new.env = new.env, proj.name = tmp_dir, xy.new.env = NULL, 
                                     selected.models = needed_predictions, compress = TRUE, 
                                     build.clamping.mask = F, do.stack = T, silent = T, 
                                     on_0_1000 = on_0_1000)
    formal_pred <- get_predictions(formal_pred, full.name = needed_predictions, 
                                   as.data.frame = ifelse(inherits(new.env, "Raster"), 
                                                          F, T))
    unlink(file.path(EM.output@sp.name, paste0("proj_", tmp_dir, 
                                              sep = "")), recursive = TRUE, force = TRUE)
  }
  ef.out <- NULL
  for (em.comp in EM.output@em.computed[which(EM.output@em.computed %in% 
                                              selected.models)]) {
    cat("\n\n\t> Projecting", em.comp, "...")
    model.tmp <- NULL
    file_name_tmp <- file.path(indiv_proj_dir, paste0(em.comp, 
                                                      output.format))
    BLM=BIOMOD_LoadModels(EM.output)
    step1=EM.output@em.models
    step2=step1[[BLM]]@model
    BIOMOD_LoadModels(EM.output, full.name = em.comp, as = "model.tmp")
    if (inherits(formal_pred, "Raster")) {
      ef.tmp <- predict(model.tmp, formal_predictions = raster::subset(formal_pred, 
                                                                       subset = step2, drop = FALSE), on_0_1000 = on_0_1000, 
                        filename = ifelse(output.format == ".RData", 
                                          "", file_name_tmp))
    }
    else {
      cat("\n*** here!!\n")
      print(str(formal_pred))
      cat("\n***\n")
      cat(model.tmp@model)
      cat("\n*** here!!\n")
      ef.tmp <- predict(model.tmp, formal_predictions = formal_pred[, 
                                                                    model.tmp@model, drop = FALSE], on_0_1000 = on_0_1000)
    }
    if (inherits(ef.tmp, "Raster")) {
      if (do.stack) {
        if (length(ef.out)) 
          ef.out <- stack(ef.out, ef.tmp)
        else ef.out <- raster::stack(ef.tmp)
      }
      else {
        file_name_tmp <- file.path(indiv_proj_dir, paste0(em.comp, 
                                                         output.format, sep = ""))
        if (output.format == ".RData") {
          save(ef.tmp, file = file_name_tmp, compress = compress)
        }
        saved.files <- c(saved.files, file_name_tmp)
      }
    }
    else {
      ef.out <- cbind(ef.out, ef.tmp)
    }
  }
  proj_out@models.projected <- EM.output@em.computed[which(EM.output@em.computed %in% 
                                                             selected.models)]
  if (do.stack) {
    if (inherits(ef.out, "Raster")) {
      names(ef.out) <- proj_out@models.projected
    }
    else {
      colnames(ef.out) <- proj_out@models.projected
    }
    file_name_tmp <- file.path(EM.output@sp.name, paste0("proj_", 
                                                        proj.name, sep = ""), paste0("proj_", proj.name, 
                                                                                    "_", EM.output@sp.name, "_ensemble", output.format, 
                                                                                    sep = ""))
    if (output.format == ".RData") {
      save(ef.out, file = file_name_tmp, compress = compress)
    }
    else if (inherits(ef.out, "Raster")) {
      writeRaster(ef.out, filename = file_name_tmp, overwrite = TRUE)
    }
    saved.files <- c(saved.files, file_name_tmp)
    proj_out@proj@link <- file_name_tmp
  }
  else {
    proj_out@proj@link <- saved.files
  }
  if (!is.null(ef.out)) {
    proj_out@proj@val <- ef.out
    proj_out@proj@inMemory <- TRUE
  }
  if (length(binary.meth) | length(filtered.meth)) {
    cat("\n")
    eval.meth <- unique(c(binary.meth, filtered.meth))
    thresholds <- sapply(selected.models, function(x) {
      get_evaluations(EM.output)[[x]][eval.meth, "Cutoff"]
    })
    if (!on_0_1000) {
      thresholds <- thresholds/1000
    }
    for (eval.meth in binary.meth) {
      cat("\n\t> Building", eval.meth, "binaries")
      if (!do.stack) {
        for (i in 1:length(proj_out@proj@link)) {
          file.tmp <- proj_out@proj@link[i]
          thres.tmp <- thresholds[i]
          proj_bin <- BinaryTransformation(raster(file.tmp, 
                                                  RAT = FALSE), thres.tmp)
          writeRaster(x = proj_bin, filename = sub(output.format, 
                                                   paste0("_", eval.meth, "bin", output.format, 
                                                         sep = ""), file.tmp), overwrite = TRUE)
        }
      }
      else {
        assign(x = paste0("proj_", proj.name, "_", EM.output@sp.name, 
                         "_ensemble_", eval.meth, "bin", sep = ""), 
               value = BinaryTransformation(ef.out, thresholds))
        if (output.format == ".RData") {
          save(list = paste0("proj_", proj.name, "_", 
                            EM.output@sp.name, "_ensemble_", eval.meth, 
                            "bin", sep = ""), file = file.path(EM.output@sp.name, 
                                                               paste0("proj_", proj.name, sep = ""), paste0("proj_", 
                                                                                                          proj.name, "_", EM.output@sp.name, "_ensemble_", 
                                                                                                          eval.meth, "bin", output.format, sep = "")), 
               compress = compress)
        }
        else {
          writeRaster(x = get(paste0("proj_", proj.name, 
                                    "_", EM.output@sp.name, "_ensemble_", eval.meth, 
                                    "bin", sep = "")), filename = file.path(EM.output@sp.name, 
                                                                            paste0("proj_", proj.name, sep = ""), paste0("proj_", 
                                                                                                                       proj.name, "_", EM.output@sp.name, "_ensemble_", 
                                                                                                                       eval.meth, "bin", output.format, sep = "")), 
                      overwrite = TRUE)
        }
        rm(list = paste0("proj_", proj.name, "_", EM.output@sp.name, 
                        "_ensemble_", eval.meth, "bin", sep = ""))
      }
    }
    for (eval.meth in filtered.meth) {
      cat("\n\t> Building", eval.meth, "filtered")
      if (!do.stack) {
        for (i in 1:length(proj_out@proj@link)) {
          file.tmp <- proj_out@proj@link[i]
          thres.tmp <- thresholds[i]
          filt_proj <- FilteringTransformation(raster(file.tmp, 
                                                      RAT = FALSE), thres.tmp)
          writeRaster(x = filt_proj, filename = sub(output.format, 
                                                    paste0("_", eval.meth, "filt", output.format, 
                                                          sep = ""), file.tmp), overwrite = TRUE)
        }
      }
      else {
        assign(x = paste0("proj_", proj.name, "_", EM.output@sp.name, 
                         "_ensemble_", eval.meth, "filt", sep = ""), 
               value = FilteringTransformation(ef.out, thresholds))
        if (output.format == ".RData") {
          save(list = paste0("proj_", proj.name, "_", 
                            EM.output@sp.name, "_ensemble_", eval.meth, 
                            "filt", sep = ""), file = file.path(EM.output@sp.name, 
                                                                paste0("proj_", proj.name, sep = ""), paste0("proj_", 
                                                                                                           proj.name, "_", EM.output@sp.name, "_ensemble_", 
                                                                                                           eval.meth, "filt", output.format, sep = "")), 
               compress = compress)
        }
        else {
          writeRaster(x = get(paste0("proj_", proj.name, 
                                    "_", EM.output@sp.name, "_ensemble_", eval.meth, 
                                    "filt", sep = "")), filename = file.path(EM.output@sp.name, 
                                                                             paste0("proj_", proj.name, sep = ""), paste0("proj_", 
                                                                                                                        proj.name, "_", EM.output@sp.name, "_ensemble_", 
                                                                                                                        eval.meth, "filt", output.format, sep = "")), 
                      overwrite = TRUE)
        }
        rm(list = paste0("proj_", proj.name, "_", EM.output@sp.name, 
                        "_ensemble_", eval.meth, "filt", sep = ""))
      }
    }
  }
  if (!keep.in.memory) {
    proj_out <- free(proj_out)
  }
  assign(paste0(EM.output@sp.name, ".", proj.name, ".ensemble.projection.out", 
               sep = ""), proj_out)
  save(list = paste0(EM.output@sp.name, ".", proj.name, ".ensemble.projection.out", 
                    sep = ""), file = file.path(EM.output@sp.name, paste0("proj_", 
                                                                         proj.name, sep = ""), paste0(EM.output@sp.name, ".", 
                                                                                                     proj.name, ".ensemble.projection.out", sep = "")))
  cat("\n")
  return(proj_out)
}


#####################################################################################
#####################################################################################
######################        Modelling procedure     ###############################
#####################################################################################
#####################################################################################


########## changing the name of the variables in the file created step 9 (Liste_predictors_final)
##########into the corresponding raster filenames  #######

predictors_list=read.csv2(paste0(predictor_list_path,"Liste_predictors_final.csv"),stringsAsFactors = FALSE)
for (i in 2:ncol(predictors_list)){
    for (j in 1:nrow(predictors_list)){
      if (is.na(predictors_list[j,i])==FALSE){
        predictors_list[j,i]=paste(sep="",raster_path,predictors_list$species[j],"_",as.character(predictors_list[j,i]),".tif")
        }
      }
    }
write.csv2(predictors_list,'predictors_list.csv',row.names = F)


#setting up parallel
library (doSNOW)
cl = makeCluster(10) ### here we choos to do 3 at a time but it can be change according to the computer performance
registerDoSNOW(cl)
iterations = nrow(predictors_list)
pb = txtProgressBar(max = iterations, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)
myLoadedPackages=(.packages())

#parallelizing
foreach (m = 1:nrow(predictors_list),.options.snow=opts,.packages=myLoadedPackages) %dopar% {    
  
  
setwd(parent_folder)
  
### create a directory for your temporary files, you must create one per species!!!  
tmpdir=paste0('D:/Users/190384/Documents/Model/Adrien_tmp/',m)
dir.create(tmpdir,recursive = T)
rasterOptions(tmpdir = tmpdir)
  
#section 1
##### give the name of the species to model ####
################################################
################################################

out1 = try ({
  
  myRespName1=predictors_list$species[m]  
  myRespName=str_replace_all(myRespName1," ",".")

})
if (class(out1) %in% 'try-error'){
  WD=getwd()
  setwd(sp_with_errors_folder)
  write.table(myRespName,paste0(myRespName,"_errors_section_1.txt"))
  setwd(WD)}


#section 2
#### selection of 80% of the occurrences randomly, 20% are kept for outer validation ####
########################################################################################
########################################################################################

out2 = try ({
  myRespXY=read.csv2(paste0(sep="",occ_aggregated_path,"occ_aggregated_ ",myRespName1," .csv"),stringsAsFactors = FALSE)
  evaluationXY=myRespXY[sample(nrow(myRespXY), (20/100)*nrow(myRespXY)), ] ### keep 20% of the occurence data for evaluation
  myRespXY_model=anti_join(myRespXY,evaluationXY)
  myResp=SpatialPoints(myRespXY_model) ### create a spatialPoints object with the 80% of occurrences kept for moddeling

})
if (class(out2) %in% 'try-error'){
  WD=getwd()
  setwd(sp_with_errors_folder)
  write.table(myRespName,paste0(myRespName,"_errors_section_2.txt"))
  setwd(WD)}

#### we want to model species that have more than 50 occurrences else it is ESM (cf step 11)####
####################################################################################
####################################################################################


if (nrow(myRespXY_model)>50){

  
###### loading the climate files into a rasterstack #####
##########################################################
##########################################################

out3 = try({files=NULL
for (i in 2:ncol(predictors_list)){
  if (is.na(predictors_list[m,i])==FALSE){    
    files=c(files,predictors_list[m,i])       
  }
}


myExpl <- stack(files)
})
if (class(out3) %in% 'try-error'){
  WD=getwd()
  setwd(sp_with_errors_folder)
  write.table(myRespName,paste0(myRespName,"_errors_section_3.txt"))
  setwd(WD)}


##### formating data ###############
####################################
####################################

out4 = try({myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                                 expl.var = myExpl,
                                                 resp.name = myRespName,
                                                 PA.nb.rep = 3,
                                                 PA.nb.absences = 10000,
                                                 PA.strategy = 'random')


myBiomodOption <- BIOMOD_ModelingOptions(MAXENT.Phillips=list(path_to_maxent.jar))

})
if (class(out4) %in% 'try-error'){
  WD=getwd()
  setwd(sp_with_errors_folder)
  write.table(myRespName,paste0(myRespName,"_errors_section_4.txt"))
  setwd(WD)}

#### modelling the species diustribution with standard models ####
##################################################################
##################################################################


out5 = try({setwd(parent_folder)
  myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData,
    models = c('GLM','GBM','MAXENT.Phillips'),
    models.options = myBiomodOption,
    NbRunEval=2,
    DataSplit=80,
    Prevalence=0.5,
    VarImport=3,
    models.eval.meth = c('TSS'),
    SaveObj = TRUE,
    rescal.all.models = FALSE,
    do.full.models = FALSE,
    modeling.id = paste0(myRespName,"_Standard_Modeling",sep=""))
  })
if (class(out5) %in% 'try-error'){
  WD=getwd()
  setwd(sp_with_errors_folder)
  write.table(myRespName,paste0(myRespName,"_errors_section_5.txt"))
  setwd(WD)}


#### outer validation and evaluation of both strategies to choose wich one to project ###
#########################################################################################
#########################################################################################

setwd(paste0(parent_folder,myRespName)) ### directory is species_dependant here with myRespName

### formating presence data kept for validation ####
####################################################

out10 = try({
### extracting climate values at presence points kept for validation  
evaluation_presence_variables=as.data.frame(raster::extract(myExpl,evaluationXY)) 
### putting 1 for presences in the first columns. The next columns are variable values
eval_presence_data_set=cbind('1',evaluation_presence_variables)

### selecting pseudo-absence data from the explanatory variables ####
#####################################################################

XY_for_pa=as.data.frame(randomPoints(myExpl,10000)) ### random choice of 10000 points in the climate data maps
### extracting climate values at presence points kept for validation  
evaluation_pseudo_absence_variables=as.data.frame(raster::extract(myExpl,XY_for_pa))
### putting 0 for absences in the first columns. The next columns are variable values
eval_pa_data_set=cbind('0',evaluation_pseudo_absence_variables)

names(eval_presence_data_set)[1]=myRespName
names(eval_pa_data_set)[1]=myRespName

### putt all the rows (presences and pseudo-absences) in the same data frame and keep it for later
eval_data_set=rbind(eval_presence_data_set,eval_pa_data_set)

write.csv2(eval_data_set,paste0(myRespName,"_dataset_for_outer_validation.csv"))})
if (class(out10) %in% 'try-error'){
  WD=getwd()
  setwd(sp_with_errors_folder)
  write.table(myRespName,paste0(myRespName,"_errors_section_10.txt"))
  setwd(WD)}



####### doing ensemble of standard models ###############
####### ensembling the 2*3=6 repetitions of each models ###########
####### to have a single ensemble model per type of modelling procedure (GLM, GBM, Maxent) ##############

### for GLM ### 
out6= try({
  setwd(parent_folder)
  myBiomodEMGLM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = get_built_models(myBiomodModelOut)[c(1,4,7,10,13,16)],
  em.by='all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.0),
  models.eval.meth = c('TSS'),
  prob.mean = T,
  prob.cv = F,
  prob.ci = F,
  prob.median = F,
  committee.averaging = F,
  prob.mean.weight = F,
  prob.mean.weight.decay = 'proportional',
  VarImport=3)

myLoadGLM=BIOMOD_LoadModels(myBiomodEMGLM)    
BIOMOD_LoadModels(myBiomodEMGLM,full.name=myLoadGLM,as="mod1")
Predict_GLM=predict (myExpl,mod1)
BOYCE_GLM= ecospat.boyce(Predict_GLM,evaluationXY)$Spearman.cor
  
GLMdir=paste0(parent_folder,myRespName,'/GLM_ensemble')
dir.create(GLMdir)
file.rename(paste0(parent_folder,myRespName,'/models/',myRespName,'_Standard_Modeling/',get_built_models(myBiomodEMGLM)),paste0(GLMdir,'/',get_built_models(myBiomodEMGLM)))  
file.rename(paste0(parent_folder,myRespName,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'),paste0(GLMdir,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'))  
})
if (class(out6) %in% 'try-error'){
  WD=getwd()
  setwd(sp_with_errors_folder)
  write.table(myRespName,paste0(myRespName,"_errors_section_6.txt"))
  setwd(WD)}


out7 = try({
  setwd(parent_folder)
  myBiomodEMGBM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = get_built_models(myBiomodModelOut)[c(2,5,8,11,14,17)],
  em.by='all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.0),
  models.eval.meth = c('TSS','ROC'),
  prob.mean = T,
  prob.cv = F,
  prob.ci = F,
  prob.median = F,
  committee.averaging = F,
  prob.mean.weight = F,
  prob.mean.weight.decay = 'proportional',
  VarImport=3)

myLoadGBM=BIOMOD_LoadModels(myBiomodEMGBM)    
BIOMOD_LoadModels(myBiomodEMGBM,full.name=myLoadGBM,as="mod2")
Predict_GBM=predict (myExpl,mod2)
BOYCE_GBM= ecospat.boyce(Predict_GBM,evaluationXY)$Spearman.cor  


GBMdir=paste0(parent_folder,myRespName,'/GBM_ensemble')
dir.create(GBMdir)
file.rename(paste0(parent_folder,myRespName,'/models/',myRespName,'_Standard_Modeling/',get_built_models(myBiomodEMGBM)),paste0(GBMdir,'/',get_built_models(myBiomodEMGBM)))    
file.rename(paste0(parent_folder,myRespName,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'),paste0(GBMdir,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'))  

})
if (class(out7) %in% 'try-error'){
  WD=getwd()
  setwd(sp_with_errors_folder)
  write.table(myRespName,paste0(myRespName,"_errors_section_7.txt"))
  setwd(WD)}


out8 = try({
  setwd(parent_folder)
  myBiomodEMMAXENT <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = get_built_models(myBiomodModelOut)[c(3,6,9,12,15,18)],
  em.by='all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.0),
  models.eval.meth = c('TSS'),
  prob.mean = T,
  prob.cv = F,
  prob.ci = F,
  prob.median = F,
  committee.averaging = F,
  prob.mean.weight = F,
  prob.mean.weight.decay = 'proportional',
  VarImport=3)

myLoadMAXENT=BIOMOD_LoadModels(myBiomodEMMAXENT)  
BIOMOD_LoadModels(myBiomodEMMAXENT,full.name=myLoadMAXENT,as="mod3")
Predict_MAXENT=predict (myExpl,mod3)
BOYCE_MAXENT= ecospat.boyce(Predict_MAXENT,evaluationXY)$Spearman.cor
    
MAXdir=paste0(parent_folder,myRespName,'/MAX_ensemble')
dir.create(MAXdir)  
file.rename(paste0(parent_folder,myRespName,'/models/',myRespName,'_Standard_Modeling/',get_built_models(myBiomodEMMAXENT)),paste0(MAXdir,'/',get_built_models(myBiomodEMMAXENT)))    
file.rename(paste0(parent_folder,myRespName,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'),paste0(MAXdir,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'))  

})
if (class(out8) %in% 'try-error'){
  WD=getwd()
  setwd(sp_with_errors_folder)
  write.table(myRespName,paste0(myRespName,"_errors_section_8.txt"))
  setwd(WD)}




#### Modelling the distribution with ensemble. ####
#### This time we ensemble all types of models ###
#### together with a condition on TSS value #####



TSS=as.data.frame(t(as.data.frame(get_evaluations(myBiomodModelOut))[1,c(1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61,65,69)]))
TSS_max=max(TSS$TSS)

out9= try({
  setwd(parent_folder)
  myBiomodEMENSEMBLE <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(min(0.7,TSS_max-0.01)),
  models.eval.meth = c('TSS'),
  prob.mean = F,
  prob.cv = F,
  prob.ci = F,
  prob.median = F,
  committee.averaging = F,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional',
  VarImport=3)

myLoadEMModels=BIOMOD_LoadModels(myBiomodEMENSEMBLE)
BIOMOD_LoadModels(myBiomodEMENSEMBLE,full.name=myLoadEMModels,as="mod4")
Predict_ENSEMBLE=predict (myExpl,mod4)
BOYCE_ENSEMBLE= ecospat.boyce(Predict_ENSEMBLE,evaluationXY)$Spearman.cor



ENSEMBLEdir=paste0(parent_folder,myRespName,'/Ensemble')
dir.create(ENSEMBLEdir)    
file.rename(paste0(parent_folder,myRespName,'/models/',myRespName,'_Standard_Modeling/',get_built_models(myBiomodEMENSEMBLE)),paste0(ENSEMBLEdir,'/',get_built_models(myBiomodEMENSEMBLE)))    
file.rename(paste0(parent_folder,myRespName,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'),paste0(ENSEMBLEdir,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'))  
})
if (class(out9) %in% 'try-error'){
  WD=getwd()
  setwd(sp_with_errors_folder)
  write.table(myRespName,paste0(myRespName,"_errors_section_9.txt"))
  setwd(WD)}



#### selection of the best model accordind to boyve index ###############
################################################
################################################


out15 = try({

setwd(paste0(parent_folder,myRespName))
  
BOYCE_max=max(BOYCE_GLM,BOYCE_GBM,BOYCE_MAXENT,BOYCE_ENSEMBLE)


if (BOYCE_max==BOYCE_GLM){
  myBiomodtoproject=myBiomodEMGLM
  projected_model='GLM'
  write.table(projected_model,'GLM_projected.txt')
  file.copy(paste0(GLMdir,'/',get_built_models(myBiomodEMGLM)),paste0(parent_folder,myRespName,'/models/',myRespName,'_Standard_Modeling/',get_built_models(myBiomodEMGLM)))
  file.copy(paste0(GLMdir,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'),paste0(parent_folder,myRespName,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'))
}

if (BOYCE_max==BOYCE_GBM){
  myBiomodtoproject=myBiomodEMGBM
  projected_model='GBM'
  write.table(projected_model,'GBM_projected.txt')
  file.copy(paste0(GBMdir,'/',get_built_models(myBiomodEMGBM)),paste0(parent_folder,myRespName,'/models/',myRespName,'_Standard_Modeling/',get_built_models(myBiomodEMGBM)))
  file.copy(paste0(GBMdir,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'),paste0(parent_folder,myRespName,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'))
}

if (BOYCE_max==BOYCE_MAXENT){
  myBiomodtoproject=myBiomodEMMAXENT
  projected_model='MAXENT'
  write.table(projected_model,'MAXENT_projected.txt')
  file.copy(paste0(MAXdir,'/',get_built_models(myBiomodEMMAXENT)),paste0(parent_folder,myRespName,'/models/',myRespName,'_Standard_Modeling/',get_built_models(myBiomodEMMAXENT)))
  file.copy(paste0(MAXdir,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'),paste0(parent_folder,myRespName,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'))
}


if (BOYCE_max==BOYCE_ENSEMBLE){
  myBiomodtoproject=myBiomodEMENSEMBLE
  projected_model='ENSEMBLE'
  write.table(projected_model,'ENSEMBLE_projected.txt')
  file.copy(paste0(ENSEMBLEdir,'/',get_built_models(myBiomodEMENSEMBLE)),paste0(parent_folder,myRespName,'/models/',myRespName,'_Standard_Modeling/',get_built_models(myBiomodEMENSEMBLE)))
  file.copy(paste0(ENSEMBLEdir,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'),paste0(parent_folder,myRespName,'/',myRespName,'.',myRespName,'_Standard_Modelingensemble.models.out'))
}


})
if (class(out15) %in% 'try-error'){
  WD=getwd()
  setwd(sp_with_errors_folder)
  write.table(myRespName,paste0(myRespName,"_errors_section_15.txt"))
  setwd(WD)}



#### projection ####
####################
####################

### selection of myRespName climate variable maps for projection #####
####################################################

out16 = try({env_files=list.files(path = future_climate_path ,full.names = TRUE,pattern = "\\.tif$")


sp_bioclim=NULL
for (k in 1:length(files)){
  for (l in 1:19){
    if (is.na(str_extract(files[k],paste0("_",l,".tif",sep="")))==FALSE){
      sp_bioclim=c(sp_bioclim,l)
    }
  }
}

sp_env=NULL
for (i in sp_bioclim){
  for (j in 1:527){
    if (is.na(str_extract(env_files[j],paste0("_",i,"_",sep="")))==FALSE){
      sp_env=c(sp_env,env_files[j])}  
  }
}


sp_26_C_61_80=NULL
sp_26_F_61_80=NULL
sp_26_MP_61_80=NULL
sp_26_MI_61_80=NULL
sp_45_C_61_80=NULL
sp_45_F_61_80=NULL
sp_45_MP_61_80=NULL
sp_45_MI_61_80=NULL
sp_85_C_61_80=NULL
sp_85_F_61_80=NULL
sp_85_MP_61_80=NULL
sp_85_MI_61_80=NULL
sp_present=NULL

for (file in sp_env){
  if (is.na(str_extract(file,"_aggregate5"))==FALSE){
    sp_present=c(sp_present,file)
  }
  if (is.na(str_extract(file,"2061-2080"))==FALSE){
    if (is.na(str_extract(file,"rcp26"))==FALSE){
      if (is.na(str_extract(file,"CESM1"))==FALSE){
        sp_26_C_61_80=c(sp_26_C_61_80,file)
      }
      if (is.na(str_extract(file,"MPI"))==FALSE){
        sp_26_MP_61_80=c(sp_26_MP_61_80,file)
      }
      if (is.na(str_extract(file,"MIROC5"))==FALSE){
        sp_26_MI_61_80=c(sp_26_MI_61_80,file)
      }
      if (is.na(str_extract(file,"FIO"))==FALSE){
        sp_26_F_61_80=c(sp_26_F_61_80,file)
      }
    }
    if (is.na(str_extract(file,"rcp45"))==FALSE){
      if (is.na(str_extract(file,"CESM1"))==FALSE){
        sp_45_C_61_80=c(sp_45_C_61_80,file)
      }
      if (is.na(str_extract(file,"MPI"))==FALSE){
        sp_45_MP_61_80=c(sp_45_MP_61_80,file)
      }
      if (is.na(str_extract(file,"MIROC5"))==FALSE){
        sp_45_MI_61_80=c(sp_45_MI_61_80,file)
      }
      if (is.na(str_extract(file,"FIO"))==FALSE){
        sp_45_F_61_80=c(sp_45_F_61_80,file)
      } 
    }
    if (is.na(str_extract(file,"rcp85"))==FALSE){
      if (is.na(str_extract(file,"CESM1"))==FALSE){
        sp_85_C_61_80=c(sp_85_C_61_80,file)
      }
      if (is.na(str_extract(file,"MPI"))==FALSE){
        sp_85_MP_61_80=c(sp_85_MP_61_80,file)
      }
      if (is.na(str_extract(file,"MIROC5"))==FALSE){
        sp_85_MI_61_80=c(sp_85_MI_61_80,file)
      }
      if (is.na(str_extract(file,"FIO"))==FALSE){
        sp_85_F_61_80=c(sp_85_F_61_80,file)
      }    
    }
  }
}


Projection_env_sp_26_C_61_80 <- raster::stack(sp_26_C_61_80,bands=1)
Projection_env_sp_26_F_61_80 <- raster::stack(sp_26_F_61_80,bands=1)
Projection_env_sp_26_MP_61_80 <- raster::stack(sp_26_MP_61_80,bands=1)
Projection_env_sp_26_MI_61_80 <- raster::stack(sp_26_MI_61_80,bands=1)
Projection_env_sp_45_C_61_80 <- raster::stack(sp_45_C_61_80,bands=1)
Projection_env_sp_45_F_61_80 <- raster::stack(sp_45_F_61_80,bands=1)
Projection_env_sp_45_MP_61_80 <- raster::stack(sp_45_MP_61_80,bands=1)
Projection_env_sp_45_MI_61_80 <- raster::stack(sp_45_MI_61_80,bands=1)
Projection_env_sp_85_C_61_80 <- raster::stack(sp_85_C_61_80,bands=1)
Projection_env_sp_85_F_61_80 <- raster::stack(sp_85_F_61_80,bands=1)
Projection_env_sp_85_MP_61_80 <- raster::stack(sp_85_MP_61_80,bands=1)
Projection_env_sp_85_MI_61_80 <- raster::stack(sp_85_MI_61_80,bands=1)
Projection_env_sp_present <- raster::stack(sp_present,bands=1)

Projection_env_sp_26_C_61_80@filename  <- '26_C_61_80'
Projection_env_sp_26_F_61_80@filename <- '26_F_61_80'
Projection_env_sp_26_MP_61_80@filename <- '26_MP_61_80'
Projection_env_sp_26_MI_61_80@filename <- '26_MI_61_80'
Projection_env_sp_45_C_61_80@filename <- '45_C_61_80'
Projection_env_sp_45_F_61_80@filename <- '45_F_61_80'
Projection_env_sp_45_MP_61_80@filename <- '45_MP_61_80'
Projection_env_sp_45_MI_61_80@filename <- '45_MI_61_80'
Projection_env_sp_85_C_61_80@filename <- '85_C_61_80'
Projection_env_sp_85_F_61_80@filename <- '85_F_61_80'
Projection_env_sp_85_MP_61_80@filename <- '85_MP_61_80'
Projection_env_sp_85_MI_61_80@filename <- '85_MI_61_80'
Projection_env_sp_present@filename <- 'present'



for (k in 1:length(sp_bioclim)){
  names(Projection_env_sp_present)[k]=names(myExpl)[k]
}

for (k in 1:length(sp_bioclim)){
  names(Projection_env_sp_26_C_61_80)[k]=names(myExpl)[k]
}

for (k in 1:length(sp_bioclim)){
  names(Projection_env_sp_26_F_61_80)[k]=names(myExpl)[k]
}

for (k in 1:length(sp_bioclim)){
  names(Projection_env_sp_26_MP_61_80)[k]=names(myExpl)[k]
}

for (k in 1:length(sp_bioclim)){
  names(Projection_env_sp_26_MI_61_80)[k]=names(myExpl)[k]
}

for (k in 1:length(sp_bioclim)){
  names(Projection_env_sp_45_C_61_80)[k]=names(myExpl)[k]
}

for (k in 1:length(sp_bioclim)){
  names(Projection_env_sp_45_F_61_80)[k]=names(myExpl)[k]
}

for (k in 1:length(sp_bioclim)){
  names(Projection_env_sp_45_MP_61_80)[k]=names(myExpl)[k]
}

for (k in 1:length(sp_bioclim)){
  names(Projection_env_sp_45_MI_61_80)[k]=names(myExpl)[k]
}

for (k in 1:length(sp_bioclim)){
  names(Projection_env_sp_85_C_61_80)[k]=names(myExpl)[k]
}

for (k in 1:length(sp_bioclim)){
  names(Projection_env_sp_85_F_61_80)[k]=names(myExpl)[k]
}

for (k in 1:length(sp_bioclim)){
  names(Projection_env_sp_85_MP_61_80)[k]=names(myExpl)[k]
}

for (k in 1:length(sp_bioclim)){
  names(Projection_env_sp_85_MI_61_80)[k]=names(myExpl)[k]
}



L=c(Projection_env_sp_26_C_61_80,
    Projection_env_sp_26_F_61_80,
    Projection_env_sp_26_MP_61_80,
    Projection_env_sp_26_MI_61_80,
    Projection_env_sp_45_C_61_80,
    Projection_env_sp_45_F_61_80,
    Projection_env_sp_45_MP_61_80,
    Projection_env_sp_45_MI_61_80,
    Projection_env_sp_85_C_61_80,
    Projection_env_sp_85_F_61_80,
    Projection_env_sp_85_MP_61_80,
    Projection_env_sp_85_MI_61_80,
    Projection_env_sp_present)

})
if (class(out16) %in% 'try-error'){
  WD=getwd()
  setwd(sp_with_errors_folder)
  write.table(myRespName,paste0(myRespName,"_errors_section_16.txt"))
  setwd(WD)}

########  Projection ###########
####################################

setwd(parent_folder)

out17 = try({


maskedEnvList=extrapolationMask(L[[length(L)]],
                                dat=evaluation_presence_variables,
                                rule='sd',
                                val=1,
                                printReport=FALSE)

maskedEnv=maskedEnvList$maskedEnv


### Projection for current climate

BIOMOD_EnsembleForecasting(myBiomodtoproject,
                                      new.env=maskedEnv,
                                      binary.meth = 'TSS',
                                      proj.name = paste0(myRespName,"_",L[[length(L)]]@filename,sep=''),
                                      compress = T)

Current_proj=stack(paste0(parent_folder,myRespName,'/proj_',paste0(myRespName,"_",L[[length(L)]]@filename,sep=''),'/proj_',paste0(myRespName,"_",L[[length(L)]]@filename,sep=''),'_',myRespName,'_ensemble_TSSbin.gri'))


setwd(parent_folder)

### Projection for all future climate scenarios and rcp

for (i in 1:(length(L)-1)){
  
  setwd(parent_folder)
  
  maskedEnvList=extrapolationMask(L[[i]],
                    dat=evaluation_presence_variables,
                    rule='sd',
                    val=1,
                    printReport=FALSE)
  
  maskedEnv=maskedEnvList$maskedEnv
  
  
  
BIOMOD_EnsembleForecasting(myBiomodtoproject,
                             new.env=maskedEnv,
                             binary.meth = 'TSS',
                             proj.name = paste0(myRespName,"_",L[[i]]@filename,sep=''),
                             compress = T)
  
Future_proj=stack(paste0(parent_folder,myRespName,'/proj_',paste0(myRespName,"_",L[[i]]@filename,sep=''),'/proj_',paste0(myRespName,"_",L[[i]]@filename,sep=''),'_',myRespName,'_ensemble_TSSbin.gri'))
  
Change=BIOMOD_RangeSize(Current_proj,Future_proj)
save(Change,file=paste0(parent_folder,myRespName,'/',myRespName,'_range_change_',L[[i]]@filename,'.Rdata'))

}
})
if (class(out17) %in% 'try-error'){
  WD=getwd()
  setwd(sp_with_errors_folder)
  write.table(myRespName,paste0(myRespName,"_errors_section_17.txt"))
  setwd(WD)}



} 

else {
return ('small species distribution not implemented')
}


done=myRespName
write.table(done,paste0(done,"_done.txt"))

setwd(tmpdir)
lf=list.files(tmpdir)
file.remove(lf)

gc()

}

#end parallel
registerDoSEQ()









BIOMOD_EnsembleForecasting=function (EM.output=NULL,
                                     projection.output = NULL,
                                     new.env = NULL, 
                                     xy.new.env = NULL,
                                     proj.name = NULL, 
                                     binary.meth = NULL,
                                     filtered.meth = NULL,
                                     compress = NULL) 
{
  # .bmCat("Do Ensemble Models Projections")
  # args <- list(...)
  # args_checked <- .BIOMOD_EnsembleForecasting.check.args(EM.output, 
  #                                                        projection.output, new.env, selected.models, proj.name, 
  #                                                        total.consensus = FALSE, binary.meth, filtered.meth)
  # proj.name <- args_checked$proj.name
  # selected.models <- args_checked$selected.models
  
  ## add
  selected.models = get_built_models(EM.output)
  
  
  # output.format <- args$output.format
  
  ## add
  output.format = NULL
  
  # compress <- args$compress
  # do.stack <- args$do.stack
  # keep.in.memory <- args$keep.in.memory
  # on_0_1000 <- args$on_0_1000
  
  ## add
  do.stack=NULL
  keep.in.memory=NULL
  on_0_1000=NULL
  
  
  if (is.null(output.format)) {
    if (length(projection.output)) {
      if (projection.output@type != "RasterStack") 
        output.format <- ".RData"
      else output.format <- ".grd"
    }
    else {
      if (!inherits(new.env, "Raster")) 
        output.format <- ".RData"
      else output.format <- ".grd"
    }
  }
  if (is.null(compress)){ 
    compress <- FALSE}
  if (is.null(do.stack)) {
    do.stack <- TRUE
    if (!is.null(projection.output)) {
      if (all(grepl("individual_projections", projection.output@proj@link))) {
        do.stack <- FALSE
      }
    }
  }
  if (is.null(keep.in.memory)) {
    keep.in.memory <- TRUE
    if (!is.null(projection.output)) {
      keep.in.memory <- projection.output@proj@inMemory
    }
  }
  if (is.null(xy.new.env)) {
    if (!is.null(projection.output)) {
      xy.new.env <- projection.output@xy.coord
    }
    else {
      xy.new.env <- matrix()
    }
  }
  if (is.null(on_0_1000)){
    on_0_1000 = TRUE
    # rm(list = c("args_checked", "args"))
    proj_out <- new("BIOMOD.projection.out", proj.names = proj.name, 
                    sp.name = EM.output@sp.name, expl.var.names = EM.output@expl.var.names, 
                    models.projected = selected.models, xy.coord = xy.new.env, 
                    modeling.object.id = EM.output@modeling.id)
    proj_out@modeling.object@link = EM.output@link
    proj_is_raster <- FALSE}
  if (inherits(new.env, "Raster")) {
    proj_is_raster <- TRUE
  }
  # else if (length(projection.output)) {
  #   if (inherits(projection.output@proj, "BIOMOD.stored.raster.stack")) {
  #     proj_is_raster <- TRUE
  #   }
  # }
  if (proj_is_raster) {
    proj_out@proj <- new("BIOMOD.stored.raster.stack")
  }
  # else {
  #   proj_out@proj <- new("BIOMOD.stored.array")
  #   do.stack = TRUE
  # }
  dir.create(file.path(parent_folder,EM.output@sp.name, paste0("proj_", proj.name, 
                                                                   sep = "")), showWarnings = FALSE, recursive = TRUE, 
             mode = "777")
  indiv_proj_dir <- file.path(parent_folder,EM.output@sp.name, paste0("proj_", 
                                                                          proj.name, sep = ""), "individual_projections")
  dir.create(indiv_proj_dir, showWarnings = FALSE, recursive = TRUE, 
             mode = "777")
  saved.files <- c()
  needed_predictions <- get_needed_models(EM.output, selected.models = selected.models)
  if (length(projection.output)) {
    formal_pred <- get_predictions(projection.output, full.name = needed_predictions,
                                   as.data.frame = ifelse(projection.output@type ==
                                                            "array", T, F))
  }
  else {
    tmp_dir <- paste0("Tmp", format(Sys.time(), "%s"), sep = "")
    formal_pred <- BIOMOD_Projection(modeling.output = load_stored_object(EM.output@models.out.obj), 
                                     new.env = new.env, proj.name = tmp_dir, xy.new.env = NULL, 
                                     selected.models = needed_predictions, compress = TRUE, 
                                     build.clamping.mask = F, do.stack = T, silent = T, 
                                     on_0_1000 = on_0_1000)
    formal_pred <- get_predictions(formal_pred, full.name = needed_predictions, 
                                   as.data.frame = ifelse(inherits(new.env, "Raster"), 
                                                          F, T))
    unlink(file.path(EM.output@sp.name, paste0("proj_", tmp_dir, 
                                              sep = "")), recursive = TRUE, force = TRUE)
  }
  ef.out <- NULL
  for (em.comp in EM.output@em.computed[which(EM.output@em.computed %in% 
                                              selected.models)]) {
    cat("\n\n\t> Projecting", em.comp, "...")
    model.tmp <- NULL
    file_name_tmp <- file.path(indiv_proj_dir, paste0(em.comp, 
                                                      output.format))
    BLM=BIOMOD_LoadModels(EM.output)
    step1=EM.output@em.models
    step2=step1[[BLM]]@model
    BIOMOD_LoadModels(EM.output, full.name = em.comp, as = "model.tmp")
    if (inherits(formal_pred, "Raster")) {
      ef.tmp <- predict(model.tmp, formal_predictions = raster::subset(formal_pred, 
                                                                       subset = step2, drop = FALSE), on_0_1000 = on_0_1000, 
                        filename = ifelse(output.format == ".RData", 
                                          "", file_name_tmp))
    }
    else {
      cat("\n*** here!!\n")
      print(str(formal_pred))
      cat("\n***\n")
      cat(model.tmp@model)
      cat("\n*** here!!\n")
      ef.tmp <- predict(model.tmp, formal_predictions = formal_pred[, 
                                                                    model.tmp@model, drop = FALSE], on_0_1000 = on_0_1000)
    }
    if (inherits(ef.tmp, "Raster")) {
      if (do.stack) {
        if (length(ef.out)) 
          ef.out <- stack(ef.out, ef.tmp)
        else ef.out <- raster::stack(ef.tmp)
      }
      else {
        file_name_tmp <- file.path(indiv_proj_dir, paste0(em.comp, 
                                                         output.format, sep = ""))
        if (output.format == ".RData") {
          save(ef.tmp, file = file_name_tmp, compress = compress)
        }
        saved.files <- c(saved.files, file_name_tmp)
      }
    }
    else {
      ef.out <- cbind(ef.out, ef.tmp)
    }
  }
  proj_out@models.projected <- EM.output@em.computed[which(EM.output@em.computed %in% 
                                                             selected.models)]
  if (do.stack) {
    if (inherits(ef.out, "Raster")) {
      names(ef.out) <- proj_out@models.projected
    }
    else {
      colnames(ef.out) <- proj_out@models.projected
    }
    file_name_tmp <- file.path(EM.output@sp.name, paste0("proj_", 
                                                        proj.name, sep = ""), paste0("proj_", proj.name, 
                                                                                    "_", EM.output@sp.name, "_ensemble", output.format, 
                                                                                    sep = ""))
    if (output.format == ".RData") {
      save(ef.out, file = file_name_tmp, compress = compress)
    }
    else if (inherits(ef.out, "Raster")) {
      writeRaster(ef.out, filename = file_name_tmp, overwrite = TRUE)
    }
    saved.files <- c(saved.files, file_name_tmp)
    proj_out@proj@link <- file_name_tmp
  }
  else {
    proj_out@proj@link <- saved.files
  }
  if (!is.null(ef.out)) {
    proj_out@proj@val <- ef.out
    proj_out@proj@inMemory <- TRUE
  }
  if (length(binary.meth) | length(filtered.meth)) {
    cat("\n")
    eval.meth <- unique(c(binary.meth, filtered.meth))
    thresholds <- sapply(selected.models, function(x) {
      get_evaluations(EM.output)[[x]][eval.meth, "Cutoff"]
    })
    if (!on_0_1000) {
      thresholds <- thresholds/1000
    }
    for (eval.meth in binary.meth) {
      cat("\n\t> Building", eval.meth, "binaries")
      if (!do.stack) {
        for (i in 1:length(proj_out@proj@link)) {
          file.tmp <- proj_out@proj@link[i]
          thres.tmp <- thresholds[i]
          proj_bin <- BinaryTransformation(raster(file.tmp, 
                                                  RAT = FALSE), thres.tmp)
          writeRaster(x = proj_bin, filename = sub(output.format, 
                                                   paste0("_", eval.meth, "bin", output.format, 
                                                         sep = ""), file.tmp), overwrite = TRUE)
        }
      }
      else {
        assign(x = paste0("proj_", proj.name, "_", EM.output@sp.name, 
                         "_ensemble_", eval.meth, "bin", sep = ""), 
               value = BinaryTransformation(ef.out, thresholds))
        if (output.format == ".RData") {
          save(list = paste0("proj_", proj.name, "_", 
                            EM.output@sp.name, "_ensemble_", eval.meth, 
                            "bin", sep = ""), file = file.path(EM.output@sp.name, 
                                                               paste0("proj_", proj.name, sep = ""), paste0("proj_", 
                                                                                                          proj.name, "_", EM.output@sp.name, "_ensemble_", 
                                                                                                          eval.meth, "bin", output.format, sep = "")), 
               compress = compress)
        }
        else {
          writeRaster(x = get(paste0("proj_", proj.name, 
                                    "_", EM.output@sp.name, "_ensemble_", eval.meth, 
                                    "bin", sep = "")), filename = file.path(EM.output@sp.name, 
                                                                            paste0("proj_", proj.name, sep = ""), paste0("proj_", 
                                                                                                                       proj.name, "_", EM.output@sp.name, "_ensemble_", 
                                                                                                                       eval.meth, "bin", output.format, sep = "")), 
                      overwrite = TRUE)
        }
        rm(list = paste0("proj_", proj.name, "_", EM.output@sp.name, 
                        "_ensemble_", eval.meth, "bin", sep = ""))
      }
    }
    for (eval.meth in filtered.meth) {
      cat("\n\t> Building", eval.meth, "filtered")
      if (!do.stack) {
        for (i in 1:length(proj_out@proj@link)) {
          file.tmp <- proj_out@proj@link[i]
          thres.tmp <- thresholds[i]
          filt_proj <- FilteringTransformation(raster(file.tmp, 
                                                      RAT = FALSE), thres.tmp)
          writeRaster(x = filt_proj, filename = sub(output.format, 
                                                    paste0("_", eval.meth, "filt", output.format, 
                                                          sep = ""), file.tmp), overwrite = TRUE)
        }
      }
      else {
        assign(x = paste0("proj_", proj.name, "_", EM.output@sp.name, 
                         "_ensemble_", eval.meth, "filt", sep = ""), 
               value = FilteringTransformation(ef.out, thresholds))
        if (output.format == ".RData") {
          save(list = paste0("proj_", proj.name, "_", 
                            EM.output@sp.name, "_ensemble_", eval.meth, 
                            "filt", sep = ""), file = file.path(EM.output@sp.name, 
                                                                paste0("proj_", proj.name, sep = ""), paste0("proj_", 
                                                                                                           proj.name, "_", EM.output@sp.name, "_ensemble_", 
                                                                                                           eval.meth, "filt", output.format, sep = "")), 
               compress = compress)
        }
        else {
          writeRaster(x = get(paste0("proj_", proj.name, 
                                    "_", EM.output@sp.name, "_ensemble_", eval.meth, 
                                    "filt", sep = "")), filename = file.path(EM.output@sp.name, 
                                                                             paste0("proj_", proj.name, sep = ""), paste0("proj_", 
                                                                                                                        proj.name, "_", EM.output@sp.name, "_ensemble_", 
                                                                                                                        eval.meth, "filt", output.format, sep = "")), 
                      overwrite = TRUE)
        }
        rm(list = paste0("proj_", proj.name, "_", EM.output@sp.name, 
                        "_ensemble_", eval.meth, "filt", sep = ""))
      }
    }
  }
  if (!keep.in.memory) {
    proj_out <- free(proj_out)
  }
  assign(paste0(EM.output@sp.name, ".", proj.name, ".ensemble.projection.out", 
               sep = ""), proj_out)
  save(list = paste0(EM.output@sp.name, ".", proj.name, ".ensemble.projection.out", 
                    sep = ""), file = file.path(EM.output@sp.name, paste0("proj_", 
                                                                         proj.name, sep = ""), paste0(EM.output@sp.name, ".", 
                                                                                                     proj.name, ".ensemble.projection.out", sep = "")))
  cat("\n")
  return(proj_out)
}
