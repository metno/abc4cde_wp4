#!/usr/bin/env Rscript
setwd(system("find $HOME -name download_tests.R -exec dirname {} \\;",intern=TRUE))
## To install DECM package: 
## R CMD INSTALL abc4cde_wp4/back-end 
## Requires rgdal, raster, esd (zoo, ncdf4), PCICt 
library(DECM)

opt <- list(verbose=TRUE,reference="era",it=c(1981,2010),variable="tas",
            nfiles=9,continue=FALSE,mask="coords.txt",help=FALSE)

#Function to calculate basic statistics
calculate.statistics.cmip <- function(reference="era", period=c(1981,2010), variable="tas", 
                                      nfiles=5, continue=TRUE, verbose=FALSE, mask="coords.txt") {
  shape <- get.shapefile("referenceRegions.shp")
  srex.regions <- as.character(shape$LAB)
  store.file <- paste("statistics.cmip", reference, variable, paste(period, collapse="-"), "rda", sep=".")
  store <- list()
  if(file.exists(store.file)) load(store.file)
  
  ref.file <- getReference(reference,variable)
  reference.raster <- raster(ref.file)
  
  store.name <- paste(reference,variable,sep=".")
  store[[store.name]]$spatial.sd <- c(cdo.spatSd(ref.file,period), cdo.spatSd(ref.file,period,seasonal=TRUE))
  store[[store.name]]$mean <- c(cdo.mean(ref.file,period), cdo.mean(ref.file,period,seasonal=TRUE))
  
  for(i in 1:length(srex.regions)){
    getPolCoords(i,shape=shape,destfile=mask)
    store[[ store.name ]][[ srex.regions[i] ]]$spatial.sd <- c(cdo.spatSd(ref.file,period,mask=mask), cdo.spatSd(ref.file,period,mask=mask,seasonal=TRUE))
    store[[ store.name ]][[ srex.regions[i] ]]$mean <- c(cdo.mean(ref.file,period,mask=mask), cdo.mean(ref.file,period,mask=mask,seasonal=TRUE))
  }
  
  ngcm <- length(cmip5.urls(varid=variable))
  start <- 1
  if(continue && file.exists(store.file))
    start <- as.numeric(tail(sub('.*\\.', '', names(store), perl=TRUE),n=1))+1
  if(nfiles=="all"){
    end <- ngcm
  } else {
    end <- min(start+nfiles-1,ngcm) 
  }
  
  for(i in start:end){
    X <- getGCMs(select=i,varid=variable)
    gcm.file <- X[[1]]$filename
    store.name <- paste("gcm",i,sep=".")
    store[[store.name]]$spatial.sd <- c(cdo.spatSd(gcm.file,period),cdo.spatSd(gcm.file,period,seasonal=T))
    store[[store.name]]$mean <- c(cdo.mean(gcm.file,period),cdo.mean(gcm.file,period,seasonal=T))
    store[[store.name]]$corr <- c(cdo.gridcor(gcm.file,ref.file,period),cdo.gridcor(gcm.file,ref.file,period,seasonal=T))
    for(j in 1:length(srex.regions)){
      getPolCoords(i,shape=shape,destfile=mask)
      store[[store.name]][[srex.regions[j]]]$spatial.sd <- c(cdo.spatSd(gcm.file,period,mask=mask), cdo.spatSd(gcm.file,period,mask=mask,seasonal=T))
      store[[store.name]][[srex.regions[j]]]$mean <- c(cdo.mean(gcm.file,period,mask=mask), cdo.mean(gcm.file,period,mask=mask,seasonal=T))
      store[[store.name]][[srex.regions[j]]]$corr <- c(cdo.gridcor(gcm.file,ref.file,period,mask=mask), cdo.gridcor(gcm.file,ref.file,period,mask=mask,seasonal=T))
    }
    save(file=store.file,store)
    gc()
    if(i==ngcm) return(store)
  }
  return(store)
}

calculate.statistics.cordex <- function(reference="era", period=c(1981,2010), variable="tas", 
                                        nfiles=5, continue=TRUE, verbose=FALSE, mask="coords.txt"){
  
  region <- read.csv(find.file("RegionSpecifications.csv"))
  region.id <- as.character(region$Code)
  store.file <- paste("statistics.cordex", reference, variable, paste(period, collapse="-"), "rda", sep=".")
  store <- list()
  if(file.exists(store.file)) load(store.file)
  
  ref.file <- getReference(reference,variable)
  reference.raster <- raster(ref.file)
  
  store.name <- paste(reference,variable,sep=".")
  store[[store.name]]$spatial.sd <- c(cdo.spatSd(ref.file,period), cdo.spatSd(ref.file,period,seasonal=T))
  store[[store.name]]$mean <- c(cdo.mean(ref.file,period), cdo.mean(ref.file,period,seasonal=T))
  
  for(i in 1:length(region.id)) {
    getPolCoords(i,shape=shape,destfile=mask)
    store[[ store.name ]][[ srex.regions[i] ]]$spatial.sd <- c(cdo.spatSd(ref.file,period,mask=mask), cdo.spatSd(ref.file,period,mask=mask,seasonal=T))
    store[[ store.name ]][[ srex.regions[i] ]]$mean <- c(cdo.mean(ref.file,period,mask=mask), cdo.mean(ref.file,period,mask=mask,seasonal=T))
  }
  
  ngcm <- length(cordex.urls(varid=variable))
  start <- 1
  if(continue && file.exists(store.file))
    start <- as.numeric(tail(sub('.*\\.', '', names(store), perl=T),n=1))+1
  if(nfiles=="all") {
    end <- ngcm
  } else {
    end <- min(start+nfiles-1,ngcm) 
  }
  
  for(i in start:end){
    X <- getRCMs(select=i,varid=variable)
    gcm.file <- X[[1]]$filename
    store.name <- paste("gcm",i,sep=".")
    store[[store.name]]$spatial.sd <- c(cdo.spatSd(gcm.file,period),cdo.spatSd(gcm.file,period,seasonal=T))
    store[[store.name]]$mean <- c(cdo.mean(gcm.file,period),cdo.mean(gcm.file,period,seasonal=T))
    store[[store.name]]$corr <- c(cdo.gridcor(gcm.file,ref.file,period),cdo.gridcor(gcm.file,ref.file,period,seasonal=T))
    for(j in 1:length(srex.regions)){
      getPolCoords(i,shape=shape,destfile=mask)
      store[[store.name]][[srex.regions[j]]]$spatial.sd <- c(cdo.spatSd(gcm.file,period,mask=mask), cdo.spatSd(gcm.file,period,mask=mask,seasonal=T))
      store[[store.name]][[srex.regions[j]]]$mean <- c(cdo.mean(gcm.file,period,mask=mask), cdo.mean(gcm.file,period,mask=mask,seasonal=T))
      store[[store.name]][[srex.regions[j]]]$corr <- c(cdo.gridcor(gcm.file,ref.file,period,mask=mask), cdo.gridcor(gcm.file,ref.file,period,mask=mask,seasonal=T))
    }
    save(file=store.file,store)
    gc()
    if(i==ngcm) return(store)
  }
  return(store)
}

calculate.statistics.cmip(reference=opt$reference, period=opt$it, variable=opt$variable, 
                          nfiles=opt$nfiles, continue=opt$continue, verbose=opt$verbose, 
                          mask=opt$mask)

calculate.statistics.cmip(reference=opt$reference, period=opt$it, variable="precip", 
                          nfiles=opt$nfiles, continue=opt$continue, verbose=opt$verbose, 
                          mask=opt$mask)

calculate.statistics.cmip(reference=opt$reference, period=opt$it, variable=opt$variable, 
                          nfiles=opt$nfiles, continue=opt$continue, verbose=opt$verbose, 
                          mask=opt$mask)

calculate.statistics.cmip(reference=opt$reference, period=opt$it, variable=opt$variable, 
                          nfiles=opt$nfiles, continue=opt$continue, verbose=opt$verbose, 
                          mask=opt$mask)

