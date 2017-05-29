## ABC4CDE/DEMC - R-script for prototype tool WP4
## Rasmus.Benestad@met.no  Oslo, Norway, 2017-02-14
##

require(esd) ## This code builds on the esd package: https://github.com/metno/esd
## A fancy colorscale:
#if(!require(wesanderson)) install.packages("wesanderson")
#library(wesanderson)

##  Function of an R-package that retrieves data files with CMIP5 or CORDEX results

getatt <- function(fname) {
  ## Reads and extracts the attribute information in a netCDF files and stores this in a list object## 
  ## This is part of the information stored in the metadatabase
  ncid <- nc_open(fname)
  nc_close(ncid)
  return(ncid)
}

## Generic function to retrieve climate model (CM) file from the KNMI ClimateExplorer
getCM <- function(url=NULL,destfile='CM.nc',lon=NULL,lat=NULL,force=FALSE,verbose=FALSE) {
  if(verbose) print("getCM")
  ## Retrieves the data
  if(is.null(url)) url <-
    'https://climexp.knmi.nl/CMIP5/monthly/tas/tas_Amon_ACCESS1-0_historical_000.nc'
  if (!file.exists(destfile)|force) lok <- try(download.file(url=url, destfile=destfile))
  if (inherits(loc,"try-error")) return()
  X <- retrieve(destfile,lon=lon,lat=lat,verbose=verbose)
  ## Collect information stored in the netCDF header
  cid <- getatt(destfile)
  ## Extract a time series for the area mean for 
  cid$area.mean <- aggregate.area(X,FUN='mean')
  cid$area.sd <- aggregate.area(X,FUN='sd')
  cid$url <- url
  cid$dates <- paste(range(index(X)),collapse=",")
  ## Collect information stored as model attributes
  ncid <- nc_open(destfile)
  model <- ncatt_get(ncid,0)
  nc_close(ncid)
  cid$model <- model
  cid$project_id <- cid$model$project_id
  return(cid)
}

## Specific function to retrieve GCMs
getGCMs <- function(select=1:9,varid='tas',destfile=NULL,verbose=FALSE) {
  if(verbose) print("getGCMs")
  ## Set destfile
  if(is.null(destfile)) destfile <- paste(rep('GCM',length(select)),select,'.',varid,'.nc',sep='')
  ## Get the urls
  url <- cmip5.urls(varid=varid)[select] ## Get the URLs of the 
  ## Set up a list variable to contain all the metadata in sub-lists.
  X <- list()
  for (i in seq_along(select)) {
    if(verbose) print(paste("Get gcm.",select[i],sep=''))
    X[[paste('gcm',varid,select[i],sep='.')]] <-
      getCM(url=url[i],destfile=destfile[i],verbose=verbose)
  }
  invisible(X)
}

testGCM <- function(select=1:9,varid='tas',path=NULL,verbose=FALSE) {
  if(verbose) print("testGCM")
  if(is.null(path)) path <- getwd()
  fnames <- list.files(path=path,pattern=varid,full.names = TRUE)
  X <- list()
  for (i in select) {
    if(verbose) print(fnames[i])
    x <- retrieve(fnames[i],varid=varid)
    ncid <- nc_open(fnames[i])
    nc_close(ncid)
    ncid$area.mean <- aggregate.area(x,FUN='mean')
    ncid$area.sd <- aggregate.area(x,FUN='sd')
    ncid$url <- fnames[i]
    X[[as.character(i)]] <- ncid
  }
  return(X)
}

## Specific model to retrieve RCMs
getRCMs <- function(select=1:9,varid='tas',destfile=NULL,verbose=FALSE) {
  if(verbose) print("getRCMs")
  ## Set destfiles
  if(is.null(destfile)) destfile <- paste('CM',select,'.',varid,'.nc',sep='')
  ## Get the urls
  url <- cordex.urls(varid=varid)[select]
  ## Set up a list variable to contain all the metadata
  X <- list()
  for (i in seq_along(select)) {
    if(verbose) print(paste("Get rcm.",select[i],sep=""))
    X[[paste('rcm',varid,select[i],sep='.')]] <- getCM(url=url[i],destfile=destfile[i],verbose=verbose)
  }
  invisible(X)
}

cmip5.urls <- function(experiment='rcp45',varid='tas',
                       url="http://climexp.knmi.nl/CMIP5/monthly/",#path=NULL,
		       off=FALSE,force=FALSE,verbose=FALSE) {
  if(verbose) print("cmip5.urls")
  urlfiles <- "NA"
  #if(is.null(path)) path <- getwd()
  for (iexp in experiment) {
    if(verbose) print(iexp)
    for (ivar in varid) {
      if(verbose) print(ivar)
      ## Loop on the number of experiments
      for (irun in 0:110) { ## 
        if(verbose) print(paste(irun))
        ## Update experiment number
        if (irun < 10) run.id = paste("00",as.character(irun),sep="")
        else if (irun < 100) run.id = paste("0",as.character(irun),sep="")
        else run.id <- as.character(irun)
        
        urlfile  <- paste(url,ivar,sep="")             # add var directory
        urlfile  <- paste(urlfile,ivar,sep="/")        # add v.name
        urlfile  <- paste(urlfile,"_Amon_ens_",sep="") # add text
        urlfile  <- paste(urlfile,iexp,sep="")         # add exp.name
        urlfile  <- paste(urlfile,run.id,sep="_")      # add exp ID number
        urlfile  <- paste(urlfile,".nc",sep="")        # add file ext
        urlfiles <- c(urlfiles,urlfile)
        if (verbose) print(urlfile)
      }
    } 
  }
  return(urlfiles[-1])
}

cordex.urls <- function(experiment='rcp45',varid='tas',
                        url="https://climexp.knmi.nl/CORDEX/EUR-44/mon/",#path=NULL,
			off=FALSE,force=FALSE,verbose=FALSE) {
  if(verbose) print("cordex.urls")
  urlfiles <- "NA"
  #if(is.null(path)) path <- getwd()
  for (iexp in experiment) {
    if(verbose) print(iexp)
    for (ivar in varid) {
      if(verbose) print(ivar)
      ## Loop on the number of experiments
      for (irun in 0:20) { ## 
        if(verbose) print(paste(irun))
        ## Update experiment number
        if (irun < 10) run.id = paste("00",as.character(irun),sep="")
        else if (irun < 100) run.id = paste("0",as.character(irun),sep="")
        else run.id <- as.character(irun)
        
        urlfile  <- paste(url,ivar,sep="")             # add var directory
        urlfile  <- paste(urlfile,ivar,sep="/")        # add v.name
        urlfile  <- paste(urlfile,"EUR-44_cordex",sep="_") # add text
        urlfile  <- paste(urlfile,iexp,"mon",sep="_")         # add exp.name
        urlfile  <- paste(urlfile,run.id,sep="_")      # add exp ID number
        urlfile  <- paste(urlfile,".nc",sep="")        # add file ext
        urlfiles <- c(urlfiles,urlfile)
        if (verbose) print(urlfile)
      }
    }
  }
  return(urlfiles[-1])
}
