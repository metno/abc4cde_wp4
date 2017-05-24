#!/usr/bin/env Rscript

#Set environmental variables
Sys.setenv(EXTERNAL_DATA = "/home/ubuntu/Data")
Sys.setenv(PROTOTYPE_SCRIPT = "/home/ubuntu/abc4cde_wp4/R-scripts")
Sys.setenv(PROTOTYPE_BACKEND = "/home/ubuntu/abc4cde_wp4/back-end")
Sys.setenv(PROTOTYPE_BACKEND_R = "/home/ubuntu/abc4cde_wp4/back-end/R")
Sys.setenv(PROTOTYPE_FRONTEND = "/home/ubuntu/abc4cde_wp4/front-end")
Sys.setenv(PROTOTYPE_DATA = "/home/ubuntu/abc4cde_wp4/back-end/data")

#Command-line parameters
suppressPackageStartupMessages({
  require(optparse)
  require(esd)
  require(raster)
  require(rgdal)
  require(RCurl)
  require(XML)
})

suppressPackageStartupMessages({
  source(paste(Sys.getenv("PROTOTYPE_BACKEND_R"),"cds.R",sep="/"))
})

option_list <- list(
  make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
              help = "Make the script verbose [default %default]"),
  make_option(c("-q", "--quietly"), action = "store_false",
              dest = "verbose", help = "Suppress run-time information"),
  make_option(c("-r", "--reference"), action = "store", default = "era",
              help = "Reference data set[default %default]"),
  make_option(c("-s", "--period1"), action = "store",default = 1981,
              help = "reference period start [default %default]"),
  make_option(c("-e", "--period2"), action = "store",default = 2010,
              help = "reference period end [default %default]"),
  make_option(c("-a", "--variable"), action = "store",default = "tas",
              help = "variable [default %default]"),
  make_option(c("-n", "--nfiles"), action = "store",default = 1,
              help = "number of files to be processed [default %default]"),
  make_option(c("-c", "--continue"), action = "store",default = F,
              help = "continue where left last time? [default %default]"),
  make_option(c("-m", "--mask"), action = "store",default = "coords.txt",
              help = "mask file [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

#Function to calculate basic statistics
calculate.statistics.cmip <- function(reference="era", period=c(1981,2010), variable="tas", nfiles=5,
                                  continue=T, verbose=F, mask="coords.txt"){
  
  shape <- get.shapefile("referenceRegions.shp")
  srex.regions <- as.character(shape$LAB)
  store.file <- paste("statistics.cmip", reference, variable, paste(period, collapse="-"), "Rdata", sep=".")
  store <- list()
  if(file.exists(store.file))store <- readRDS(store.file)

  ref.file <- getReference(reference,variable)
  store.name <- paste(reference,variable,sep=".")
  store[[store.name]]$spatial.sd <- c(cdo.spatSd(ref.file,period), cdo.spatSd(ref.file,period,seasonal=T))
  store[[store.name]]$mean <- c(cdo.mean(ref.file,period), cdo.mean(ref.file,period,seasonal=T))
  
  for(i in 1:length(srex.regions)){
    getPolCoords(shape,i,destfile=mask)
    store[[ store.name ]][[ srex.regions[i] ]]$spatial.sd <- c(cdo.spatSd(ref.file,period,mask=mask), cdo.spatSd(ref.file,period,mask=mask,seasonal=T))
    store[[ store.name ]][[ srex.regions[i] ]]$mean <- c(cdo.mean(ref.file,period,mask=mask), cdo.mean(ref.file,period,mask=mask,seasonal=T))
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
    gcm.file <- get.name(i,variable)
    if(!file.exists(gcm.file)) download.file(cmip5.urls(i,variable), destfile=paste(system("echo $PROTOTYPE_DATA"),gcm.file,sep="/"))
    store.name <- paste("gcm",i,sep=".")
    store[[store.name]]$spatial.sd <- c(cdo.spatSd(gcm.file,period),cdo.spatSd(gcm.file,period,seasonal=T))
    store[[store.name]]$mean <- c(cdo.mean(gcm.file,period),cdo.mean(gcm.file,period,seasonal=T))
    store[[store.name]]$corr <- c(cdo.gridcor(gcm.file,ref.file,period),cdo.gridcor(gcm.file,ref.file,period,seasonal=T))
    for(j in 1:length(srex.regions)){
      getPolCoords(shape,i,destfile=mask)
      store[[store.name]][[srex.regions[j]]]$spatial.sd <- c(cdo.spatSd(gcm.file,period,mask=mask), cdo.spatSd(gcm.file,period,mask=mask,seasonal=T))
      store[[store.name]][[srex.regions[j]]]$mean <- c(cdo.mean(gcm.file,period,mask=mask), cdo.mean(gcm.file,period,mask=mask,seasonal=T))
      store[[store.name]][[srex.regions[j]]]$corr <- c(cdo.gridcor(gcm.file,ref.file,period,mask=mask), cdo.gridcor(gcm.file,ref.file,period,mask=mask,seasonal=T))
    }
    saveRDS(store,store.file)
    gc()
    if(i==ngcm)return(invisible(store))
  }
  return(invisible(store))
}

store <- calculate.statistics.cmip(opt$reference,c(opt$period1,opt$period2),opt$variable,opt$nfiles,
                          opt$continue,opt$verbose,opt$mask)
