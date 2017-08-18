#!/usr/bin/env Rscript
setwd(system("find $HOME -name download_tests.R -exec dirname {} \\;",intern=TRUE))
## To install DECM package: 
## R CMD INSTALL abc4cde_wp4/back-end 
## Requires rgdal, raster, esd (zoo, ncdf4), PCICt 
library(DECM)

#Function to calculate basic statistics
calculate.statistics.cmip <- function(reference="era", period=c(1981,2010), variable="tas", 
                                      nfiles=5, continue=TRUE, verbose=FALSE, mask="coords.txt") {
  shape <- get.shapefile("referenceRegions.shp")
  srex.regions <- as.character(shape$LAB)
  if(max(period)>2015) reference <- NULL
  if(!is.null(reference)) {
    ref.file <- getReference(reference,variable)
    if(!is.character(ref.file)) {
      reference <- NULL
      print("Warning! Reference file not found. Continuing without reference data.")
    }
  }
  if(!is.null(reference)) {
    store.file <- paste("statistics.cmip", reference, variable, paste(period, collapse="-"), "rda", sep=".")
  } else {
    store.file <- paste("statistics.cmip", variable, paste(period, collapse="-"), "rda", sep=".")
  }
  store <- list()
  if(file.exists(store.file)) load(store.file)
  
  units <- NULL
  if(!is.null(reference)) {
    reference.raster <- raster(ref.file)
    store.name <- paste(reference,variable,sep=".")
    store[[store.name]]$spatial.sd <- c(cdo.spatSd(ref.file,period), cdo.spatSd(ref.file,period,monthly=TRUE))
    store[[store.name]]$mean <- c(cdo.mean(ref.file,period), cdo.mean(ref.file,period,monthly=TRUE))
    if(variable=="pr") {
      if(max(abs(store[[store.name]]$mean),na.rm=TRUE)<0.001) {
        units <- "kg m-2 s-1"
      } else {
        units <- "mm/day"
      }
    } else if (variable=="tas") {
      if(max(abs(store[[store.name]]$mean),na.rm=TRUE)>273) {
        units <- "K"
      } else {
        units <- "degrees~Celsius"
      }
    }
    
    for(i in 1:length(srex.regions)){
      getPolCoords(i,shape=shape,destfile=mask)
      store[[ store.name ]][[ srex.regions[i] ]]$spatial.sd <- c(cdo.spatSd(ref.file,period,mask=mask), cdo.spatSd(ref.file,period,mask=mask,monthly=TRUE))
      store[[ store.name ]][[ srex.regions[i] ]]$mean <- c(cdo.mean(ref.file,period,mask=mask), cdo.mean(ref.file,period,mask=mask,monthly=TRUE))
    }
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
    store[[store.name]]$spatial.sd <- c(cdo.spatSd(gcm.file,period),cdo.spatSd(gcm.file,period,monthly=T))
    store[[store.name]]$mean <- c(cdo.mean(gcm.file,period),cdo.mean(gcm.file,period,monthly=T))
    if(!is.null(reference)) store[[store.name]]$corr <- c(cdo.gridcor(gcm.file,ref.file,period),cdo.gridcor(gcm.file,ref.file,period,monthly=T))
    for(j in 1:length(srex.regions)){
      getPolCoords(j,shape=shape,destfile=mask)
      store[[store.name]][[srex.regions[j]]]$spatial.sd <- c(cdo.spatSd(gcm.file,period,mask=mask), cdo.spatSd(gcm.file,period,mask=mask,monthly=T))
      store[[store.name]][[srex.regions[j]]]$mean <- c(cdo.mean(gcm.file,period,mask=mask), cdo.mean(gcm.file,period,mask=mask,monthly=T))
      if(!is.null(reference)) store[[store.name]][[srex.regions[j]]]$corr <- c(cdo.gridcor(gcm.file,ref.file,period,mask=mask), cdo.gridcor(gcm.file,ref.file,period,mask=mask,monthly=T))
    }
    save(file=store.file,store)
    gc()
    if(i==ngcm) return(store)
  }
  attr(store,"variable") <- variable
  if(variable=="pr") {
    if(max(abs(store[[1]]$mean),na.rm=TRUE)<0.001) {
      attr(store,"unit") <- "kg m-2 s-1"
    } else {
      attr(store,"unit") <- "mm/day"
    }
  } else if (variable=="tas") {
    if(max(abs(store[[1]]$mean),na.rm=TRUE)>273) {
      attr(store,"unit") <- "K"
    } else {
      attr(store,"unit") <- "degrees~Celsius"
    }
  }
  save(file=store.file,store)
  return(store)
}

calculate.statistics.cordex <- function(reference="era", period=c(1981,2010), variable="tas", 
                                        nfiles=5, continue=TRUE, verbose=FALSE, mask="coords.txt"){
  
  region <- read.csv(find.file("RegionSpecifications.csv"))
  region.id <- as.character(region$Code)
  if(max(period)>2015) reference <- NULL
  if(!is.null(reference)) {
    store.file <- paste("statistics.cordex", reference, variable, paste(period, collapse="-"), "rda", sep=".")
  } else {
    store.file <- paste("statistics.cordex", variable, paste(period, collapse="-"), "rda", sep=".")
  }
  
  store <- list()
  if(file.exists(store.file)) load(store.file)
  
  if(!is.null(reference)) {
    ref.file <- getReference(reference,variable)
    reference.raster <- raster(ref.file)
  
    store.name <- paste(reference,variable,sep=".")
    store[[store.name]]$spatial.sd <- c(cdo.spatSd(ref.file,period), cdo.spatSd(ref.file,period,monthly=T))
    store[[store.name]]$mean <- c(cdo.mean(ref.file,period), cdo.mean(ref.file,period,monthly=T))
  
    #for(i in 1:length(region.id)) {
    #  getPolCoords(i,shape=shape,destfile=mask)
    #  store[[ store.name ]][[ srex.regions[i] ]]$spatial.sd <- c(cdo.spatSd(ref.file,period,mask=mask), cdo.spatSd(ref.file,period,mask=mask,monthly=T))
    #  store[[ store.name ]][[ srex.regions[i] ]]$mean <- c(cdo.mean(ref.file,period,mask=mask), cdo.mean(ref.file,period,mask=mask,monthly=T))
    #}
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
    store[[store.name]]$spatial.sd <- c(cdo.spatSd(gcm.file,period),cdo.spatSd(gcm.file,period,monthly=T))
    store[[store.name]]$mean <- c(cdo.mean(gcm.file,period),cdo.mean(gcm.file,period,monthly=T))
    if(!is.null(reference)) store[[store.name]]$corr <- c(cdo.gridcor(gcm.file,ref.file,period),cdo.gridcor(gcm.file,ref.file,period,monthly=T))
    #for(j in 1:length(srex.regions)){
    #  getPolCoords(j,shape=shape,destfile=mask)
    #  store[[store.name]][[srex.regions[j]]]$spatial.sd <- c(cdo.spatSd(gcm.file,period,mask=mask), cdo.spatSd(gcm.file,period,mask=mask,monthly=T))
    #  store[[store.name]][[srex.regions[j]]]$mean <- c(cdo.mean(gcm.file,period,mask=mask), cdo.mean(gcm.file,period,mask=mask,monthly=T))
    #  if(!is.null(reference)) store[[store.name]][[srex.regions[j]]]$corr <- c(cdo.gridcor(gcm.file,ref.file,period,mask=mask), cdo.gridcor(gcm.file,ref.file,period,mask=mask,monthly=T))
    #}
    save(file=store.file,store)
    gc()
    if(i==ngcm) return(store)
  }
  return(store)
}

opt <- list(verbose=TRUE,reference="era",it=c(1981,2010),variable="tas",
            nfiles=9,continue=FALSE,mask="coords.txt",help=FALSE)
for (varid in c("tas","pr")) {
  calculate.statistics.cmip(reference=opt$reference, period=opt$it, variable=varid, 
                            nfiles=opt$nfiles, continue=opt$continue, verbose=opt$verbose, 
                            mask=opt$mask)
  for (it in list(c(2021,2050),c(2071,2100))) {
    calculate.statistics.cmip(reference=NULL, period=it, variable=varid, 
                              nfiles=opt$nfiles, continue=opt$continue, verbose=opt$verbose, 
                              mask=opt$mask)
  }
}

#calculate.statistics.cordex(reference=opt$reference, period=opt$it, variable=opt$variable, 
#                          nfiles=opt$nfiles, continue=opt$continue, verbose=opt$verbose, 
#                          mask=opt$mask)

