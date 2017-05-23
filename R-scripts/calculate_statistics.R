#!/usr/bin/env Rscript

#Command-line parameters
suppressPackageStartupMessages({
  require(optparse)
})


setwd(system("find $HOME -name test_bed.R -exec dirname {} \\;",intern=T))
source("../back-end/R/cds.R")
require(esd)
require(raster)
require(rgdal)

store.statistics.cmip(reference="era", period=c(1981,2010), variable="tas", nfiles="all", continue=T, verbose=F)

#Function to calculate spatial statistics for Taylor diagram
store.statistics.cmip <- function(reference="era", period=c(1981,2010), variable="tas", nfiles=5, continue=T, verbose=F, mask="coords.txt"){
  shape <- get.shapefile("referenceRegions.shp")
  srex.regions <- as.character(shape$LAB)
  store.file <- paste("statistics.cmip", reference, variable, paste(period, collapse="-"), "Rdata", sep=".")
  store <- list()
  if(file.exists(store.file))store <- readRDS(store.file)

  ref.file <- getReference(reference,variable)
  reference.raster <- raster(ref.file)
  
  store.name <- paste(reference,variable,sep=".")
  store[[store.name]]$spatial.sd <- c(getSpatialSd(ref.file,period), getSpatialSd(ref.file,period,seasonal=T))
  store[[store.name]]$mean <- c(getMean(ref.file,period), getMean(ref.file,period,seasonal=T))
  
  for(i in 1:length(srex.regions)){
    getPolCoords(shape,i,destfile=mask)
    store[[ store.name ]][[ srex.regions[i] ]]$spatial.sd <- c(getSpatialSd(ref.file,period,mask=mask), getSpatialSd(ref.file,period,mask=mask,seasonal=T))
    store[[ store.name ]][[ srex.regions[i] ]]$mean <- c(getMean(ref.file,period,mask=mask), getMean(ref.file,period,mask=mask,seasonal=T))
  }
  
  ngcm <- length(cmip5.urls(varid=variable))
  start <- 1
  if(continue && file.exists(store.file))
    start <- as.numeric(tail(sub('.*\\.', '', names(store), perl=T),n=1))+1
  if(nfiles=="all"){
    end <- ngcm
  }else{
    end <- min(start+nfiles-1,ngcm) 
  }
  
  for(i in start:end){
    gcm.file <- getGCM(i,variable)
    if(!file.exists(gcm.file)) download.file(cmip5.urls(i,variable), destfile=paste(system("echo $PROTOTYPE_DATA"),gcm.file,sep="/"))
    store.name <- paste("gcm",i,sep=".")
    store[[store.name]]$spatial.sd <- c(getSpatialSd(gcm.file,period),getSpatialSd(gcm.file,period,seasonal=T))
    store[[store.name]]$mean <- c(getMean(gcm.file,period),getMean(gcm.file,period,seasonal=T))
    store[[store.name]]$corr <- c(getPatternCor(gcm.file,ref.file,period),getPatternCor(gcm.file,ref.file,period,seasonal=T))
    for(j in 1:length(srex.regions)){
      getPolCoords(shape,i,destfile=mask)
      store[[store.name]][[srex.regions[j]]]$spatial.sd <- c(getSpatialSd(gcm.file,period,mask=mask), getSpatialSd(gcm.file,period,mask=mask,seasonal=T))
      store[[store.name]][[srex.regions[j]]]$mean <- c(getMean(gcm.file,period,mask=mask), getMean(gcm.file,period,mask=mask,seasonal=T))
      store[[store.name]][[srex.regions[j]]]$corr <- c(getPatternCor(gcm.file,ref.file,period,mask=mask), getPatternCor(gcm.file,ref.file,period,mask=mask,seasonal=T))
    }
    saveRDS(store,store.file)
    gc()
    if(i==ngcm)return(store)
  }
  return(store)
}

#---------------------------
#GetReference
getReference <- function(reference,variable){
  path <- system("echo $PROTOTYPE_DATA",intern=T)
  file.name <- switch(paste(reference,variable,sep="."),
         era.tas="era-interim_monthly_1979-2016_tas.2.5deg.nc",
         era.pr="era-interim_monthly_1979-2016_pr.2.5deg.nc",
         cfsr.tas="cfsr_tmp2m_mon.nc",
         cfsr.pr="cfsr_prate_mon.nc")
  invisible(paste(path,file.name,sep="/"))
}

#---------------------------
#GetGCM
getGCM <- function(number,variable){
  path <- system("echo $PROTOTYPE_DATA",intern=T)
  file.name <- paste(paste("GCM",number,sep=""),variable,"nc",sep=".")
  invisible(paste(path,file.name,sep="/"))
}

#---------------------------
#getPolCoords
getPolCoords <- function(shape,region,destfile=mask){
  if(is.character(region))region <- which(as.character(shape$LAB)==region)
  pol.coords <- coordinates(shape@polygons[[region]]@Polygons[[1]])
  write(t(pol.coords),file=destfile,ncolumns = 2)
}

#---------------------------
#GetMean
getMean <- function(model,period=c(1981,2010),mask=NULL,seasonal=F){
  
  commands <- c("-fldmean","-timmean","-selyear")
  input <- c("","",paste(period,collapse="/"))
  
  if(!is.null(mask)){
    commands <- append(commands,"-maskregion",after=2)
    input <- append(input,mask,after=2) 
  }
  if(seasonal){
    commands <- replace(commands,commands=="-timmean","-yseasmean")
  }

  out.file <- "tmp.nc"
  cdo.command(commands,input,model,out.file)
  
  command <- ("output")
  input <- c("")
  
  out <- as.numeric(cdo.command(command,input,out.file,NULL,intern=T))
  if(seasonal){
    names(out) <- c("djf","mam","jja","son")
  }else{
    names(out) <- "ann"
  } 
  
  if(out>200)out <- out-273.15
  system(paste("rm",out.file,sep=" "))
  invisible(out)
}

#---------------------------
#GetTimeSd
getTimeSd <- function(model,period=c(1981,2010),mask=NULL,seasonal=F){
  
  commands <- c("-timstd","-fldmean","-ymean","-selyear")
  input <- c("","",paste(period,collapse="/"))
    
  if(!is.null(mask)){
    commands <- append(commands,"-maskregion",after=3)
    input <- append(input,mask,after=3) 
  }
  if(seasonal){
    commands <- replace(commands,commands=="-ymean","-seasmean")
  }

  out.file <- "tmp.nc"
  cdo.command(commands,input,model,out.file)
  
  command <- ("output")
  input <- c("")
  
  out <- as.numeric(cdo.command(command,input,out.file,NULL,intern=T))
  if(seasonal){
    names(out) <- c("djf","mam","jja","son")
  }else{
    names(out) <- "ann"
  }
  system(paste("rm",out.file,sep=" "))
  invisible(out)
}
  
getSpatialSd <- function(model,period=c(1981,2010),mask=NULL,seasonal=F){

  commands <- c("-fldstd","-timmean","-selyear")
  input <- c("","",paste(period,collapse="/"))
    
  if(!is.null(mask)){
    commands <- append(commands,"-maskregion",after=2)
    input <- append(input,mask,after=2) 
  }
  if(seasonal){
    commands <- replace(commands,commands=="-timmean","-yseasmean")
  }
  
  out.file <- "tmp.nc"
  cdo.command(commands,input,model,out.file)
  
  command <- ("output")
  input <- c("")
  out <- as.numeric(cdo.command(command,input,out.file,NULL,intern=T))
  if(seasonal){
    names(out) <- c("djf","mam","jja","son")
  }else{
    names(out) <- "ann"
  }
  system("rm tmp.nc")
  invisible(out)
}

#---------------------------
#getPatternCor
getPatternCor <- function(model.file,reference.file,period=c(1981,2010),mask=NULL,seasonal=F){
  
  commands <- c("-timavg","-selyear")
  input <- c("",paste(period,collapse="/"))
  
  if(!is.null(mask)){
    commands <- append(commands,"-maskregion",after=2)
    input <- append(input,mask,after=2) 
  }
  if(seasonal){
    commands <- replace(commands,commands=="-timavg","-yseasavg")
  }
  
  out.file <- "tmp.nc"
  cdo.command(commands,input,model.file,out.file)
  
  out.file <- "tmp2.nc"
  cdo.command(commands,input,reference.file,out.file)
  
  commands <- c("fldcor")
  input <- c("")
  in.file <- c("tmp.nc tmp2.nc")
  out.file <- "tmp_cor.nc"
  cdo.command(commands,input,in.file,out.file)
  
  command <- ("output")
  input <- c("")
  out <- as.numeric(cdo.command(command,input,out.file,NULL,intern=T))
  if(seasonal){
    names(out) <- c("djf","mam","jja","son")
  }else{
    names(out) <- "ann"
  }
  system("rm tmp.nc tmp2.nc tmp_cor.nc")
  invisible(out)
}
