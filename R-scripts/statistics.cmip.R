#!/usr/bin/env Rscript
#Script that calculates statistics for a given period.
#
#Olle Räty 2017

get_scriptpath <- function() {
  # location of script can depend on how it was invoked:
  # source() and knit() put it in sys.calls()
  path <- NULL
  if(Sys.getenv("RSTUDIO") == "1"){
    path <- dirname(rstudioapi::getActiveDocumentContext()$path)
  }else if(!is.null(sys.calls())) {
    # get name of script - hope this is consisitent!
    path <- as.character(sys.call(1))[2] 
    # make sure we got a file that ends in .R, .Rmd or .Rnw
    if (grepl("..+\\.[R|Rmd|Rnw]", path, perl=TRUE, ignore.case = TRUE) )  {
      return(path)
    } else { 
      message("Obtained value for path does not end with .R, .Rmd or .Rnw: ", path)
    }
  } else{
    # Rscript and R -f put it in commandArgs
    args <- commandArgs(trailingOnly = FALSE)
  }
  return(path)
}

sourcePrototypeFile <- function(file,current.path){
  print(file)
  repository.root <- gsub("(abc4cde_wp4).*","\\1",current.path)
  print(paste("find",repository.root,"-name",file,sep=" "))
  full.path <-system(paste("find",repository.root,"-name",file,sep=" "),intern=T)
  print(full.path)
  source(full.path)
}

#Set environmental variables
current_path <- get_scriptpath()
setwd(current_path)
sourcePrototypeFile("cds.R",current_path)
suppressPackageStartupMessages({
  require(optparse)
  require(esd)
})

Sys.setenv(EXTERNAL_DATA="/home/ubuntu/Data/Data")

#Command-line parameters
option_list <- list(
  make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
              help = "Make the script verbose [default %default]"),
  make_option(c("-q", "--quietly"), action = "store_false",
              dest = "verbose", help = "Suppress run-time information"),
  make_option(c("-r", "--reference"), action = "store", default = "era",
              help = "Reference data set. If set to NULL, statistics for GCM simulations are calculated only[default %default]"),
  make_option(c("-s", "--period1"), action = "store",default = 2071,
              help = "First yeart of the calculation period [default %default]"),
  make_option(c("-e", "--period2"), action = "store",default = 2100,
              help = "Last year of the calculation period [default %default]"),
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
EXTERNAL_DATA = "/home/ubuntu/Data"

#Function to calculate basic statistics
calculate.statistics.cmip <- function(reference="era", period=c(1981,2010), variable="tas", nfiles=5,
                                  continue=T, verbose=F, mask="coords.txt"){
  
  if(is.null(reference)){
    store.file <- paste("statistics.cmip", reference, variable, paste(period, collapse="-"), "Rdata", sep=".") 
  } else{
    store.file <- paste("statistics.cmip", variable, paste(period, collapse="-"), "Rdata", sep=".")
  }
  shape <- get.shapefile("referenceRegions.shp")
  srex.regions <- as.character(shape$LAB)
  store <- list()
  if(file.exists(store.file))store <- readRDS(store.file)

  if(!is.null(reference)){
    ref.file <- getReference(reference,variable,path=Sys.getenv("EXTERNAL_DATA"))
    store.name <- paste(reference,variable,sep=".")
    store[[store.name]]$spatial.sd <- c(cdo.spatSd(ref.file,period), cdo.spatSd(ref.file,period,seasonal=T))
    store[[store.name]]$mean <- c(cdo.mean(ref.file,period), cdo.mean(ref.file,period,seasonal=T))
    
    for(i in 1:length(srex.regions)){
      getPolCoords(shape,i,destfile=mask)
      store[[ store.name ]][[ srex.regions[i] ]]$spatial.sd <- c(cdo.spatSd(ref.file,period,mask=mask), cdo.spatSd(ref.file,period,mask=mask,seasonal=T))
      store[[ store.name ]][[ srex.regions[i] ]]$mean <- c(cdo.mean(ref.file,period,mask=mask), cdo.mean(ref.file,period,mask=mask,seasonal=T))
    } 
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
    if(!file.exists(gcm.file)) download.file(cmip5.urls(i,variable), destfile=paste(Sys.getenv("EXTERNAL_DATA"),gcm.file,sep="/"))
    store.name <- paste("gcm",i,sep=".")
    store[[store.name]]$spatial.sd <- c(cdo.spatSd(gcm.file,period),cdo.spatSd(gcm.file,period,seasonal=T))
    store[[store.name]]$mean <- c(cdo.mean(gcm.file,period),cdo.mean(gcm.file,period,seasonal=T))
    if(!is.null(reference))store[[store.name]]$corr <- c(cdo.gridcor(gcm.file,ref.file,period),cdo.gridcor(gcm.file,ref.file,period,seasonal=T))
    for(j in 1:length(srex.regions)){
      getPolCoords(shape,j,destfile=mask)
      store[[store.name]][[srex.regions[j]]]$spatial.sd <- c(cdo.spatSd(gcm.file,period,mask=mask), cdo.spatSd(gcm.file,period,mask=mask,seasonal=T))
      store[[store.name]][[srex.regions[j]]]$mean <- c(cdo.mean(gcm.file,period,mask=mask), cdo.mean(gcm.file,period,mask=mask,seasonal=T))
      if(!is.null(reference))store[[store.name]][[srex.regions[j]]]$corr <- c(cdo.gridcor(gcm.file,ref.file,period,mask=mask), cdo.gridcor(gcm.file,ref.file,period,mask=mask,seasonal=T))
    }
    saveRDS(store,store.file)
  }
  return(invisible(store))
}

store <- calculate.statistics.cmip(NULL,c(opt$period1,opt$period2),opt$variable,opt$nfiles,
                          opt$continue,opt$verbose,opt$mask)
