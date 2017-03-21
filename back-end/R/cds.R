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
  ## KMP 2017-03-13: Not all important information is stored in the netCDF header.
  ##                 Collect info about model (some can also be found in object X)
  ncid <- nc_open(destfile)
  model <- ncatt_get(ncid,0)
  nc_close(ncid)
  ## Extract a time series for the area mean for 
  cid$area.mean <- aggregate.area(X,FUN='mean')
  cid$url <- url
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
    rcmnames <- c(rcmnames,attr(paste('appendix.',i,sep=''),'model_id'))
    #paste('appendix.',i,sep='') <- NULL
    clim[[paste('rcm.',i+1,sep='')]] <- map.field(attr(X,paste('appendix.',i,sep='')))
  }
  attr(Z,'mean') <- clim
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
  class(Z) <- c('dsensemble','eof','zoo')
  ceof <- Z
  save(ceof,file=paste('ceof.rcm.',it,'.rda',sep=''))
}



