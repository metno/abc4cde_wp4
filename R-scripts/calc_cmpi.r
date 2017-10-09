library(DECM)
source("/home/ubuntu/abc4cde_wp4/back-end/R/get.R") #Did not manage to install the latest version of DECM!
source("/home/ubuntu/abc4cde_wp4/back-end/R/regions.R") #Did not manage to install the latest version of DECM!

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

#Calculate weights for the weighted average of rms, assuming monthly time step.
calc.mon.weights <- function(lon,lat){
  
  weights <- array(NA,dim=c(12,length(lon),length(lat)))
  time.weights <- c(31,28,31,30,31,30,31,31,30,31,30,31)/365
  lat.weights <- rep(cos(pi*lat/180),length(lon))
  dim(lat.weights) <- c(length(lat),length(lon))
  lat.weights <- t(lat.weights)
  image(lat.weights)
  for(i in 1:length(lat)){
    weights[,,i] <- time.weights
  }
  
  for(i in 1:12){
    weights[i,,] <- weights[i,,]*lat.weights
  }
  
  return(weights)
}

calculate.cmpi <- function(reference="era",period=c(1981,2010), variable="tas", nfiles=4,
                                      continue=F, verbose=F){

  shape <-  get.shapefile("../back-end/data/SREX_regions/referenceRegions.shp",with.path = T)
  external <- "/home/ubuntu/Data/Data/" #COMPUTER SPECIFIC, CHANGE!!!!!
    
  store <- list()
  store.file <- paste("statistics.cmip", reference, variable, paste(period, collapse="-"), "rda", sep=".")
  if(file.exists(store.file)) load(store.file)
  
  #Pre-process ERA-interim, if necessary
  era.file <- paste(external,"era-interim_monthly_1979-2017_",variable,".2.5deg.nc",sep="")
  era.mulc <- paste(paste(external,"era.mulc",sep=""),variable,"nc",sep=".")
  era.mon.file <- paste(paste(external,"era.monmean",sep=""),variable,"nc",sep=".")
  
  if(variable=="pr"){
    if(!file.exists(era.mulc)) cdo.command("mulc",1000,era.file,era.mulc)
  }else{
    if(!file.exists(era.mulc)) era.mulc <- era.file
  }
  
  if(!file.exists(era.mon.file)) cdo.command(c("-ymonmean","-selyear"),c("",paste(period,collapse="/")),
             infile=era.mulc,outfile=era.mon.file)
  era <- coredata(retrieve(era.mon.file))
  
  #Calculate weights only once
  r <- raster(era.mon.file)
  lon <- xFromCol(r)
  lat <- yFromRow(r)
  weights <- calc.mon.weights(lon,lat)
  
  #Check which files are processed
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
    store.name <- paste("gcm",i,sep=".")
    gcm.file <- paste(paste(external,"GCM",i,sep=""),variable,"nc",sep=".")
    if(!file.exists(gcm.file))getGCMs(i,varid=variable,destfile=paste(external,gcm.file,sep=""))
    gcm.mon.file <- paste(external,"gcm.monmean.nc",sep="")
    cdo.command(c("-ymonmean","-selyear"),c("",paste(period,collapse="/")),
                infile=gcm.file,outfile=gcm.mon.file)
    gcm <- coredata(retrieve(gcm.mon.file))
    dim(gcm) <- dim(era) <- c(12,length(longitude(gcm)),length(latitude(gcm)))
    store[[store.name]]$rms <- sqrt(sum(weights*(gcm-era)^2)/sum(weights))
    for(region in levels(shape$LAB)){
      polygon <- shape[which(levels(shape$LAB)==region),]
      mask <- gen.mask.srex(destfile=gcm.file,mask.polygon=polygon)
      dim(gcm) <- dim(era) <- c(12,length(longitude(era))*length(latitude(era)))
      gcm.masked <- mask.zoo(gcm,mask)
      era.masked <- mask.zoo(era,mask)
      dim(gcm.masked) <- dim(era.masked) <- c(12,length(longitude(gcm)),length(latitude(gcm)))
      store[[store.name]][[region]]$rms <- sqrt(sum(weights*(gcm.masked-era.masked)^2,na.rm=T)/sum(weights[!is.na(gcm.masked)]))
    }
    file.remove(gcm.mon.file)
  }
  
  median.rms <- median(unlist(lapply(store, "[","rms")))  
  for(i in start:end){
    store[[store.name]]$e <- (store[[store.name]]$rms-median.rms)/median.rms
    for(region in levels(shape$LAB)){
      median.rms <- median(unlist(lapply(lapply(store, "[",region),function(x) x[[region]]$rms)))
      store[[store.name]][[region]]$e <- (store[[store.name]][[region]]$rms-median.rms)/median.rms
    }
  }
  save(file=store.file,store)
  file.remove(era.mon.file)
}

period <- c(1981,2010)
setwd(get_scriptpath())
variable <- "tas"
nfiles <- 4
calculate.cmpi(period=c(1981,2010), variable="tas", nfiles=4,
                           continue=F, verbose=F)
