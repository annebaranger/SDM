library(dismo)
library(ecospat)
library(raster)
library(stringr)
library(rJava)

####### package to install at first time #####
##############################################
##############################################

Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre1.8.0_241')
Sys.setenv(PATH='C:/Program Files/Java/jre1.8.0_241/bin/')

library(devtools)
install_github("johnbaums/rmaxent")


############### used fonction  #####################
####################################################
####################################################

#### project the model (use a maxent object as argument)
project <- function(lambdas, newdata, return_lfx=FALSE, mask, quiet=FALSE) {
  if(!missing(mask)) {
    if(!methods::is(mask, 'RasterLayer')) {
      stop('mask should be a RasterLayer object')
    } else {
      if(!methods::is(newdata, 'Raster')) {
        stop('If mask is provided, newdata should be a Raster object with the',
             'same dimensions, extent, and CRS.')
      }
      if(!raster::compareRaster(mask, newdata, stopiffalse=FALSE))
        stop('If mask is provided, newdata should be a Raster object with the',
             'same dimensions, extent, and CRS.')
    }
  }
  if(!is.lambdas(lambdas)) lambdas <- parse_lambdas(lambdas)
  meta <- lambdas[-1]
  lambdas <- lambdas[[1]]
  is_cat <- unique(
    gsub('\\(|==.*\\)', '', lambdas[lambdas$type=='categorical', 'feature']))
  nms <- unique(unlist(strsplit(lambdas$var[lambdas$lambda != 0], ',')))
  clamp_limits <- data.table::data.table(lambdas[lambdas$type=='linear', ])
  lambdas <- lambdas[lambdas$lambda != 0, ]
  lambdas <- split(lambdas, c('other', 'hinge')[
    grepl('hinge', lambdas$type) + 1])
  if (is(newdata, 'RasterStack') | is(newdata, 'RasterBrick') | is(newdata, 'Raster')) {
    pred_raw <- pred_logistic <- pred_cloglog <- pred_lfx <- raster::raster(newdata)
    if(!missing(mask)) {
      newdata <- raster::mask(newdata, mask)
    }
    newdata <- data.table::as.data.table(newdata[])
  }
  if (is.matrix(newdata)) newdata <- data.table::as.data.table(newdata)
  if (is.list(newdata) & !is.data.frame(newdata)) {
    if(length(unique(sapply(newdata, length))) != 1) 
      stop('newdata was provided as a list, but its elements vary in length.')
    newdata <- data.table::as.data.table(newdata)
  }
  if (!is.data.frame(newdata))
    stop('newdata must be a list, data.table, data.frame, matrix, RasterStack,',
         'or RasterBrick.')
  if (!data.table::is.data.table(newdata)) 
    newdata <- data.table::as.data.table(newdata)
  if(!all(nms %in% names(newdata))) {
    stop(sprintf('Variables missing in newdata: %s', 
                 paste(setdiff(nms, colnames(newdata)), collapse=', ')))
  }
  if (any(!names(newdata) %in% nms)) {
    newdata <- newdata[, setdiff(names(newdata), nms) := NULL]
  }
  na <- !complete.cases(newdata)
  newdata <- newdata[!na]
  
  # Clamp newdata
  invisible(lapply(setdiff(names(newdata), is_cat), function(x) {
    clamp_max <- clamp_limits[clamp_limits$feature==x, max]
    clamp_min <- clamp_limits[clamp_limits$feature==x, min]
    newdata[, c(x) := pmax(pmin(get(x), clamp_max), clamp_min)]
  }))
  
  k_hinge <- if('hinge' %in% names(lambdas)) nrow(lambdas$hinge) else 0
  k_other <- if('other' %in% names(lambdas)) nrow(lambdas$other) else 0
  k <- k_hinge + k_other
  
  txt <- sprintf('\rCalculating contribution of feature %%%1$dd of %%%1$dd', 
                 nchar(k))
  lfx <- numeric(nrow(newdata))
  lfx_all <- setNames(vector('list', sum(sapply(lambdas, nrow))),
                      unlist(lapply(lambdas[2:1], function(x) x$feature)))
  
  if(k_other > 0) {
    for (i in seq_len(k_other)) {
      if(!quiet) cat(sprintf(txt, i, k))
      x <- with(newdata, eval(parse(text=lambdas$other$feature[i])))
      # clamp feature
      x <- pmin(pmax(x, lambdas$other$min[i]), lambdas$other$max[i])
      x01 <- (x - lambdas$other$min[i]) / 
        (lambdas$other$max[i] - lambdas$other$min[i])
      lfx_all[[i]] <- lambdas$other$lambda[i] * x01
      lfx <- lfx + lfx_all[[i]]
    }
    rm(x, x01)
  }
  
  if(k_hinge > 0) {
    for (i in seq_len(nrow(lambdas$hinge))) {
      if(!quiet) cat(sprintf(txt, k_other + i, k))
      x <- with(newdata, get(sub("'|`", "", lambdas$hinge$feature[i])))
      x01 <- (x - lambdas$hinge$min[i]) / (lambdas$hinge$max[i] - lambdas$hinge$min[i])
      if (lambdas$hinge$type[i]=='reverse_hinge') {
        lfx_all[[k_other + i]] <- 
          lambdas$hinge$lambda[i] * (x < lambdas$hinge$max[i]) * (1-x01)
      } else {
        lfx_all[[k_other + i]] <- 
          lambdas$hinge$lambda[i] * (x >= lambdas$hinge$min[i]) * x01
      }
      lfx <- lfx + lfx_all[[k_other + i]]
    }
    rm(x, x01)
  }
  
  ln_raw <- lfx - meta$linearPredictorNormalizer - log(meta$densityNormalizer)
  raw <- exp(ln_raw)
  logit <- meta$entropy + ln_raw
  cloglog <- 1 - exp(-exp(meta$entropy) * raw)
  logistic <- stats::plogis(logit)
  
  #linpred <- rep(NA_real_, length(na))
  #linpred[!na] <- lfx
  if(exists('pred_raw', inherits=FALSE)) {
    pred_raw[which(!na)] <- raw
    pred_logistic[which(!na)] <- logistic
    pred_cloglog[which(!na)] <- cloglog
    pred_lfx[which(!na)] <- lfx
    out <- list(prediction_raw=pred_raw,
                prediction_logistic=pred_logistic,
                prediction_cloglog=pred_cloglog)
    if(isTRUE(return_lfx)) {
      lfx_each <- lapply(lfx_all, function(x) {
        r <- raster(pred_raw)
        r[which(!na)] <- x
        r
      })
      out <- c(out, 
               list(prediction_lfx=pred_lfx,
                    lfx_all=lfx_each))
    } 
  } else {
    prediction_raw <- prediction_logistic <- prediction_cloglog <-
      rep(NA_real_, length(na))
    prediction_raw[!na] <- raw
    prediction_logistic[!na] <- logistic
    prediction_cloglog[!na] <- cloglog
    #prediction_lfx[!na] <- lfx
    out <- list(prediction_raw=prediction_raw,
                prediction_logistic=prediction_logistic,
                prediction_cloglog=prediction_cloglog)
  }
  return(out)
}

### associated function for maxent and lambdas objects
is.lambdas <- function(x) is(x, 'lambdas')


parse_lambdas <- function(lambdas) {
  if(methods::is(lambdas, 'MaxEnt')) {
    lambdas <- lambdas@lambdas
  } else {
    lambdas <- readLines(lambdas)
  }
  con <- textConnection(lambdas)
  n <- utils::count.fields(con, ',', quote='')
  close(con)
  meta <- stats::setNames(lapply(strsplit(lambdas[n==2], ', '), 
                                 function(x) as.numeric(x[2])),
                          sapply(strsplit(lambdas[n==2], ', '), '[[', 1))
  lambdas <- stats::setNames(data.frame(do.call(
    rbind, strsplit(lambdas[n==4], ', ')), stringsAsFactors=FALSE),
    c('feature', 'lambda', 'min', 'max'))
  lambdas[, -1] <- lapply(lambdas[, -1], as.numeric)
  lambdas$feature <- sub('=', '==', lambdas$feature)
  lambdas$feature <- sub('<', '<=', lambdas$feature)
  lambdas$type <- factor(sapply(lambdas$feature, function(x) {
    switch(gsub("\\w|\\.|-|\\(|\\)", "", x),
           "==" = 'categorical',
           "<=" = "threshold",
           "^" = "quadratic",
           "*" = "product", 
           "`" = "reverse_hinge",
           "'" = 'forward_hinge',
           'linear')
  }))
  vars <- gsub("\\^2|\\(.*<=|\\((.*)==.*|`|\\'|\\)", "\\1", lambdas$feature)
  lambdas$var <- sub('\\*', ',', vars)
  l <- c(list(lambdas=lambdas[, c(1, 6, 2:5)]), meta)
  class(l) <- 'lambdas'
  l
}


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
    masks=stack(lapply(vars,function(x){
      e.min=min(dat[,x],na.rm=T)
      e.max=max(dat[,x],na.rm=T)
      e.sd=sd(dat[,x],na.rm=T)
      e.min.lim=e.min-val*e.sd
      e.max.lim=e.max+val*e.sd
      mask1=thisEnv[[x]] >= e.min.lim & thisEnv[[x]] <= e.max.lim
      values(mask1)[values(mask1)==0]=NA
      mask1
    }))
    names(masks)=vars
    #masks=raster::trim(masks)
    maskedEnv=stack(lapply(vars,function(x){
      raster::mask(thisEnv[[x]],masks[[x]],maskValue=0)
    }))
    names(maskedEnv)=vars
  }
  
  # make a report of which layers lead to masking
  if(printReport){
    nonNAPreMask=apply(values(thisEnv),2,function(x) sum(!is.na(x)))
    nonNAMasked=apply(values(maskedEnv),2,function(x) sum(!is.na(x)))
    print(data.frame(nonNAPreMask=nonNAPreMask,nonNAMasked=nonNAMasked, numCellsMasked=nonNAPreMask-nonNAMasked))
  }
  return(list(masks=masks,maskedEnv=maskedEnv))
  #})
  #return(out)
}



### where have you saved your predictors_list obtained step 9
predictor_list_path="D:/Users/190384/Documents/Model/"


### where have you saved your climate rasters obtained step 7
raster_path="D:/Write/Noirault/Data/climate_extraction_per_spe/"

### where are your aggregated occurences data obtained step 8
occ_aggregated_path="D:/Write/Noirault/Data/occ_aggregated/"

### parent folder for model results per species
parent_folder='D:/Users/190384/Documents/Model/modele/'

### where are your future climate rasters
future_climate_path="D:/Write/Noirault/Data/climat_tif/"

### species with errors folder
sp_with_errors_folder="D:/Users/190384/Documents/Model/modele/Species_with_errors/"

########## changing the name of the variables in the file created step 9 (Liste_predictors_final)
##########into the corresponding raster filenames  #######

predictors_list=read.csv2(paste0(predictor_list_path,"Liste_predictors_final.csv"),stringsAsFactors = FALSE)

for (i in 2:ncol(predictors_list)){
  for (j in 1:nrow(predictors_list)){
    if (is.na(predictors_list[j,i])==FALSE){
      predictors_list[j,i]=paste0(sep="",raster_path,predictors_list$species[j],"_",as.character(predictors_list[j,i]),".tif")
    }
  }
}

###############################################
#### Modelling procedure for small species ####
###############################################
#### Small species = If there are less than 50 occurences remaining for modelling #####
#### i.e if 80% of the records represents less than 50 #####
##### (cause 20% are kept for evaluation) ################
##########################################################


for (j in 1:nrow(predictors_list)){


myRespName=predictors_list$species[j]  
myRespXY=read.csv2(paste0(sep="",occ_aggregated_path,"occ_aggregated_ ",myRespName," .csv"),stringsAsFactors = FALSE)


out1 = try ({


if (nrow(myRespXY)*0.8<51) {

  ### create a output directory for each species
  output.dir = paste0(parent_folder,myRespName)
  dir.create(output.dir)
  setwd(output.dir)  

#### create a raster stack with all the climate variable rasters ####
files=NULL
for (i in 2:ncol(predictors_list)){
  if (is.na(predictors_list[j,i])==FALSE){    
    files=c(files,predictors_list[j,i])       
  }
}
myExpl <- stack(files) 

#### extracting climate variable values at occurence points
presence_data_climate=as.data.frame(raster::extract(myExpl,myRespXY))
presence_data=cbind(cbind('1',presence_data_climate))

#### selecting random pseudo-absences and ###
### extracting climate variable values at these points ####
XY_for_pa=as.data.frame(randomPoints(myExpl, 10000))
pseudo_absence_data_climate=as.data.frame(raster::extract(myExpl,XY_for_pa))
pa_data=cbind('0',pseudo_absence_data_climate)

names(presence_data)[1]='pa'
names(pa_data)[1]='pa'

### putting PA and presence data in the same data frame
data_set=rbind(presence_data,pa_data)

### env = table with climate values at each point
env  <- data_set[,2:ncol(data_set)]

# number of absences
sel.pa <- data_set$pa==0
sum(sel.pa)

# number of presences
sel.occ <- data_set$pa==1
sum(sel.occ)

### env.occ = env data where presence records
### env.bg = env data where PA where chosen
env.occ <- env[sel.occ,]
env.bg  <- env[sel.pa ,]

resp <- c(rep(1,nrow(env.occ)),rep(0,nrow(env.bg)))

# =================================================================================
# 10 runs using a random 80:20 split of data for training and testing ESM vs Maxent
# =================================================================================

eval.test.st  <-  NULL    # for evaluation results of Maxent standard
eval.test.esm <-  NULL    # for evaluation results of Ensemble of Small Model


for (i in 1:5) {  
  
  # split data for training and testing (80:20)
  mykfold <- kfold(resp, k=5, by=resp)
  sel.train  <- mykfold != 1
  sel.test   <- mykfold == 1
  resp.train <- resp[sel.train]
  resp.test  <- resp[sel.test]
  
  # standard Maxent
  st <- maxent(x=env[sel.train,], p=resp[sel.train], args=c("-P","nohinge", "nothreshold", "noproduct"))
  pred.st <- predict(st,env)    
  auc.test <- evaluate(p=pred.st[sel.test][resp.test==1], a=pred.st[sel.test ][resp.test ==0])@auc
  auc.test.st <- evaluate(p=pred.st[sel.test][resp.test==1], a=pred.st[sel.test ][resp.test ==0])
  boyce.test <- ecospat.boyce(na.omit(pred.st[sel.test]), pred.st[sel.test][resp.test==1], PEplot=F)$Spearman.cor
  eval.test.st <- rbind(eval.test.st, cbind(auc=auc.test,boyce=boyce.test))
  
  # build Ensemble of Small Model as weighted average of bivariate Maxent models
  weight.vec <- vector()                 # vector with weights of single bivariate models
  pred.biva.sum <- rep(0, length(resp))  # for ensemble prediction of bivariate models 
  
  for(m in 1:(dim(env)[2]-1)){           # build all bivariate predictor combinations
    for(n in (m+1):dim(env)[2]){
      
      # calibrate bivariate maxent model
      me.biva <- maxent(x=env[sel.train,c(m,n)], p=resp.train, args=c("-P","nohinge", "nothreshold", "noproduct"))  
      
      # prediction
      pred.biva <- predict(me.biva, env)
      
      # evaluate bivariate model and calculate its weight for ensemble prediction
      eval.train.biva <- evaluate(p=pred.biva[sel.train][resp.train==1], a=pred.biva[sel.train][resp.train==0])
      weight <- eval.train.biva@auc*2-1  # use Somers' D as weighting factor
      if(weight < 0) {weight <- 0}
      weight.vec <- c(weight.vec, weight)
      
      # build weighted sum of bivariate predictions
      pred.biva.sum <- pred.biva.sum + (pred.biva * weight)
      
    }
  }

save.image(file=paste0(myRespName,'_modelling.RData'))
    
  # calculate ESM prediction
  pred.esm <- pred.biva.sum/sum(weight.vec)
  
  # evaluate ESM
  auc.test <- evaluate(p=pred.esm[sel.test ][resp.test ==1], a=pred.esm[sel.test ][resp.test ==0])@auc
  auc.test.esm <- evaluate(p=pred.esm[sel.test ][resp.test ==1], a=pred.esm[sel.test ][resp.test ==0])
  boyce.test <- ecospat.boyce(pred.esm[sel.test],pred.esm[sel.test][resp.test==1], PEplot=F)$Spearman.cor     
  eval.test.esm <- rbind(eval.test.esm , cbind(auc=auc.test,boyce=boyce.test))
  
  print(i)

}


save.image(file=paste0(myRespName,'_modelling.RData'))


### save results of model evaluation into tables 
eval.test.ESM=as.data.frame(eval.test.esm)
eval.test.ST=as.data.frame(eval.test.st)
write.table(eval.test.ESM,paste0(myRespName,'_AUC_BOYCE_for_ESM.csv'),sep=";",row.names = F)
write.table(eval.test.ESM,paste0(myRespName,'_AUC_BOYCE_for_standard_Maxent.csv'),sep=";",row.names = F)

#### selection of the best model to project it to future climate

BOYCE_ESM=mean(eval.test.ESM[,2])
BOYCE_ST=mean(eval.test.ST[,2])

if (BOYCE_ST<BOYCE_ESM){
  myModeltoProject=me.biva
  save(myModeltoProject,file='projected.model_esm.RData')
  auc=auc.test.esm
} else {
  myModeltoProject=st
  save(myModeltoProject,file='projected.model_st.RData')
  auc=auc.test.st}


##### getting future climate maps for the selected variables (different scenarios and RCP)
env_files=list.files(path = future_climate_path, full.names = TRUE,pattern = "\\.tif$")


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


Projection_env_sp_26_C_61_80 <- stack(sp_26_C_61_80,bands=1)
Projection_env_sp_26_F_61_80 <- stack(sp_26_F_61_80,bands=1)
Projection_env_sp_26_MP_61_80 <- stack(sp_26_MP_61_80,bands=1)
Projection_env_sp_26_MI_61_80 <- stack(sp_26_MI_61_80,bands=1)
Projection_env_sp_45_C_61_80 <- stack(sp_45_C_61_80,bands=1)
Projection_env_sp_45_F_61_80 <- stack(sp_45_F_61_80,bands=1)
Projection_env_sp_45_MP_61_80 <- stack(sp_45_MP_61_80,bands=1)
Projection_env_sp_45_MI_61_80 <- stack(sp_45_MI_61_80,bands=1)
Projection_env_sp_85_C_61_80 <- stack(sp_85_C_61_80,bands=1)
Projection_env_sp_85_F_61_80 <- stack(sp_85_F_61_80,bands=1)
Projection_env_sp_85_MP_61_80 <- stack(sp_85_MP_61_80,bands=1)
Projection_env_sp_85_MI_61_80 <- stack(sp_85_MI_61_80,bands=1)
Projection_env_sp_present <- stack(sp_present,bands=1)


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


i=1

#### PROJECTIONS for present and future climate

for (i in 1:length(L)){
  
  
  maskedEnvList=extrapolationMask(L[[i]],
                                  dat=env.bg,
                                  rule='sd',
                                  val=1,
                                  printReport=FALSE)
  
  maskedEnv=maskedEnvList$maskedEnv
  

  x=project(myModeltoProject,newdata = maskedEnv ,return_lfx = F,quiet = F)
  

  ##### give a understandable name to the file and save it as a raster
  nRname = names (x)
  lapply (nRname , function (n){
    print (n)
    writeRaster(x[[n]],
                filename = paste0(output.dir,'/',myRespName,'_',n,'_',L[[i]]@filename,'.tif'),
                overwrite=T, datatype=  'FLT4S',options="COMPRESS=LZW")
  })
  
##### Choice of a threshold to make binary projections ####
  threshold=threshold(auc,'equal_sens_spec')
 
### if the projected presence probability is higher than the threshold, put a 1
###  else put a 0
  pres=c(threshold,1,1)
  abs=c(0,threshold,0)
  rcl=as.data.frame(rbind(pres,abs))
  Binary_proj=reclassify(x[["prediction_logistic"]],rcl)
  
  writeRaster(Binary_proj, filename = paste0(output.dir,'/',myRespName,'_',n,'_',L[[i]]@filename,'binary_projection.tif'),
              overwrite=T, datatype=  'FLT4S',options="COMPRESS=LZW")
  
}
}

  })
if (class(out1) %in% 'try-error'){
  WD=getwd()
  setwd(sp_with_errors_folder)
  write.table(myRespName,paste0(myRespName,"_errors_section_1.txt"))
  setwd(WD)}
}
