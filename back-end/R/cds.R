## ABC4CDE/DEMC - R-script for prototype tool WP4
## Rasmus.Benestad@met.no  Oslo, Norway, 2017-02-14
##

require(esd) ## This code builds on the esd package: https://github.com/metno/esd

##  Function of an R-package that retrieves data files with CMIP or CORDEX results

getatt <- function(fname) {
  ## Reads and extracts the attribute information in a netCDF files and stores this in a list object## 
  ## This is part of the information stored in the metadatabase
  ncid <- nc_open(fname)
  nc_close(ncid)
  return(ncid)
}

## Generic function to retrieve climate model (CM) file from the KNMI ClimateExplorer
getCM <- function(url='https://climexp.knmi.nl/CMIP5/monthly/tas/tas_Amon_ACCESS1-0_historical_000.nc',
                  destfile='CM.nc',lon=NULL,lat=NULL,force=FALSE) {
  ## Retrieves the data 
  if (!file.exists(destfile)|force) lok <- try(download.file(url=url,destfile = destfile)) 
  if (inherits(loc,"try-error")) return()
  X <- retrieve(destfile,lon=lon,lat=lat)
  ## Collect information stored in the netCDF header
  cid <- getatt(destfile)
  ## Extract a time series for the area mean for 
  cid$area.mean <- aggregate.area(X,FUN='mean')
  cid$url <- url
  cid$dates <- paste(range(index(X)),collapse=",")
  ## KMP 2017-03-13: Not all important information is stored in the netCDF header.
  ##                 Collect info about model (some can also be found in object X)
  ncid <- nc_open(destfile)
  model <- ncatt_get(ncid,0)
  nc_close(ncid)
  cid$model <- model
  return(cid)
}

## Specific function to retrieve GCMs
getGCMs <- function(select=1:3,varid='tas',destfile=NULL) {
  ## Set destfile
  if(is.null(destfile)) destfile <- paste(rep('GCM',length(select)),select,'.nc',sep='')
  ## Get the urls
  url <- cmip5.urls(varid=varid)[select] ## Get the URLs of the 
  ## Set up a list variable to contain all the metadata in sub-lists.
  X <- list()
  for (i in seq_along(select)) X[[paste('gcm',select[i],sep='.')]] <- getCM(url=url[i],destfile=destfile[i])
  return(X)
}

testGCM <- function(select=1:3,varid='tas',path='~/storeB/CMIP5.monthly/rcp45/') {
  fnames <- list.files(path=path,pattern=varid,full.names = TRUE)
  X <- list()
  for (i in select) {
    print(fnames[i])
    x <- retrieve(fnames[i],varid=varid)
    ncid <- nc_open(fnames[i])
    nc_close(ncid)
    ncid$area.mean <- aggregate.area(x,FUN='mean')
    ncid$url <- fnames[i]
    X[[as.character(i)]] <- ncid
  }
  return(X)
}


## Specific model to retrieve RCMs
getRCMs <- function(select=1:3,varid='tas',destfile=NULL) {
  ## Set destfiles
  if(is.null(destfile)) destfile <- paste(rep('CM',length(select)),select,'.nc',sep='')
  ## Get the urls
  ## KMP 2017-03-21: added function cordex.urls in cmip5.download.R
  url <- cordex.urls()[select]
  #url <- paste('https://climexp.knmi.nl/CORDEX/EUR-44/mon/',varid,'/',varid,'_EUR-44_cordex_rcp45_mon_00',select,'.nc',sep='')
  #url <- sub('0000','000',url)
  ## Set up a list variable to contain all the metadata
  X <- list()
  for (i in seq_along(select)) X[[paste('rcm',select[i],sep='.')]] <- getCM(url=url[i],destfile=destfile[i])
  return(X)
}

## Compute the common EOFs for GCMs and save the results for the front-end
commonEOFS.gcm <- function(select=1:3,varid='tas',destfile=NULL,
                           it='annual',is=NULL) {
  if(is.null(destfile)) destfile <- paste(rep('GCM',length(select)),select,'.nc',sep='')
  getGCMs(select=select,varid=varid,destfile=destfile)
  X <- NULL
  for (fname in destfile) {
    x <- retrieve(fname)
    if (!is.null(it)) {
      if (tolower(it)=='annual') x <- annual(x) else
                                 x <- subset(x,it=it,is=is)
    }
    if (is.null(X)) X <- x else X <- combine(X,x)
  }
  ceof <- EOF(X)
  plot(ceof)
  ## Need to reformat the ceof-object to fit the set-up for the R-shiny app in the front-end.
  ## Need to reformat the ceof-object to fit the set-up for the R-shiny app in the front-end.
  x1 <- coredata(ceof); attributes(x1) <- NULL; dim(x1) <- dim(ceof)
  Z <- list(info='CORDEX runs',eof=ceof,rcm.1=zoo(x1,order.by=index(ceof)))
  clim <- list(rcm.1=map.field(X,plot=FALSE))
  rcmnames <- attr(ceof,'model_id')
  for (i in 1:attr(ceof,'n.apps')) {
    x1 <- coredata(attr(ceof,paste('appendix.',i,sep='')))
    #attributes(x1) <- NULL; dim(x1) <- dim(attr(ceof,paste('appendix.',i,sep='')))
    Z[[paste('rcm.',i+1,sep='')]] <- zoo(x1,order.by=index(paste('appendix.',i,sep='')))
    rcmnames <- c(rcmnames,attr(x1,'model_id'))
    #paste('appendix.',i,sep='') <- NULL
    clim[[paste('rcm.',i+1,sep='')]] <- map.field(attr(X,paste('appendix.',i,sep='')))
  }
  attr(Z,'mean') <- clim
  attr(Z,'model_id') <- rcmnames
  class(Z) <- c('dsensemble','eof','zoo')
  ceof <- Z
  save(ceof,file=paste('ceof.gcm.',it,'.rda',sep=''))
}

## Compute the common EOFs for RCMs save the results for the front-end
commonEOFS.rcm <- function(select=1:3,varid='tas',destfile=NULL,
                           it='annual',is=NULL) {
  if(is.null(destfile)) destfile <- paste(rep('CM',length(select)),select,'.nc',sep='')
  getRCMs(select=select,varid=varid,destfile=destfile)
  X <- NULL
  for (fname in destfile) {
    x <- retrieve(fname)
    if (!is.null(it)) {
      if (tolower(it)=='annual') x <- annual(x) else
        x <- subset(x,it=it,is=is)
    }
    if (is.null(X)) X <- x else X <- combine(X,x)
    
  }
  ceof <- EOF(X)
  plot(ceof)
  
  ## Need to reformat the ceof-object to fit the set-up for the R-shiny app in the front-end.
  x1 <- coredata(ceof); attributes(x1) <- NULL; dim(x1) <- dim(ceof)
  Z <- list(info='CORDEX runs',eof=ceof,rcm.1=zoo(x1,order.by=index(ceof)))
  clim <- list(rcm.1=map.field(X,plot=FALSE))
  rcmnames <- attr(ceof,'model_id')
  for (i in 1:attr(ceof,'n.apps')) {
    x1 <- coredata(attr(ceof,paste('appendix.',i,sep='')))
    #attributes(x1) <- NULL; dim(x1) <- dim(attr(ceof,paste('appendix.',i,sep='')))
    Z[[paste('rcm.',i+1,sep='')]] <- zoo(x1,order.by=index(paste('appendix.',i,sep='')))
    rcmnames <- c(rcmnames,attr(paste('appendix.',i,sep=''),'model_id'))
    #paste('appendix.',i,sep='') <- NULL
    clim[[paste('rcm.',i+1,sep='')]] <- map.field(attr(X,paste('appendix.',i,sep='')))
  }
  attr(Z,'mean') <- clim
  attr(Z,'model_id') <- rcmnames
  class(Z) <- c('dsensemble','eof','zoo')
  ceof <- Z
  save(ceof,file=paste('ceof.rcm.',it,'.rda',sep=''))
}


cmip5.urls <- function(experiment='rcp45',varid='tas',
                       url="http://climexp.knmi.nl/CMIP5/monthly/", 
                       path=NULL,off=FALSE,force=FALSE,verbose=FALSE) {
  urlfiles <- "NA"
  if(verbose) print("cmip5.urls")
  if(is.null(path)) path <- getwd()
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
        
        ## Create output directory for the climate experiment
        path.exp <- file.path(path,experiment[grep(iexp,experiment)],
                              fsep = .Platform$file.sep)
        if (!file.exists(path.exp)) dir.create(path.exp)
        if (verbose) print(path.exp[grep(iexp,experiment)])
        ## Define the output file
        destfile <- paste(path.exp,varid[grep(ivar,varid)],sep="/") 
        destfile <- paste(destfile,"_Amon_ens_",sep="")
        destfile <- paste(destfile,iexp,sep="")
        destfile <- paste(destfile,run.id,sep="_")
        destfile <- paste(destfile,".nc",sep="")
        
        if (!file.exists(destfile) | force) {
          ## Update output filename with attributes:
          urlfile  <- paste(url,ivar,sep="")             # add var directory
          urlfile  <- paste(urlfile,ivar,sep="/")        # add v.name
          urlfile  <- paste(urlfile,"_Amon_ens_",sep="") # add text
          urlfile  <- paste(urlfile,iexp,sep="")         # add exp.name
          urlfile  <- paste(urlfile,run.id,sep="_")      # add exp ID number
          urlfile  <- paste(urlfile,".nc",sep="")        # add file ext
        }
        urlfiles <- c(urlfiles,urlfile)
        if (verbose) print(urlfile)
      }
      
      
    } # End for   
  }
  return(urlfiles[-1])
}



cordex.urls <- function(experiment='rcp45',varid='tas',
                        url="https://climexp.knmi.nl/CORDEX/EUR-44/mon/", 
                        path=NULL,off=FALSE,force=FALSE,verbose=FALSE) {
  urlfiles <- "NA"
  if(verbose) print("cordex.urls")
  if(is.null(path)) path <- getwd()
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
        
        ## Create output directory for the climate experiment
        path.exp <- file.path(path,experiment[grep(iexp,experiment)],
                              fsep = .Platform$file.sep)
        if (!file.exists(path.exp)) dir.create(path.exp)
        if (verbose) print(path.exp[grep(iexp,experiment)])
        ## Define the output file
        destfile <- paste(path.exp,varid[grep(ivar,varid)],sep="/") 
        destfile <- paste(destfile,"EUR-44_cordex",sep="_")
        destfile <- paste(destfile,iexp,"mon",sep="_")
        destfile <- paste(destfile,run.id,sep="_")
        destfile <- paste(destfile,".nc",sep="")
        
        if (!file.exists(destfile) | force) {
          ## Update output filename with attributes:
          urlfile  <- paste(url,ivar,sep="")             # add var directory
          urlfile  <- paste(urlfile,ivar,sep="/")        # add v.name
          urlfile  <- paste(urlfile,"EUR-44_cordex",sep="_") # add text
          urlfile  <- paste(urlfile,iexp,"mon",sep="_")         # add exp.name
          urlfile  <- paste(urlfile,run.id,sep="_")      # add exp ID number
          urlfile  <- paste(urlfile,".nc",sep="")        # add file ext
        }
        urlfiles <- c(urlfiles,urlfile)
        if (verbose) print(urlfile)
      }
      
    } # End for   
  }
  return(urlfiles[-1])
}

