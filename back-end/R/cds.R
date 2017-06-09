## ABC4CDE/DEMC - R-script for prototype tool WP4
## Rasmus.Benestad@met.no  Oslo, Norway, 2017-02-14
##


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
  if (file.exists(destfile) & !force) {
    X <- try(retrieve(destfile,lon=lon,lat=lat,verbose=verbose))
    if (inherits(X,"try-error")) force <- TRUE # If downloaded file is incomplete, force new download
  }
  if (!file.exists(destfile) | force) {
    lok <- try(download.file(url=url, destfile=destfile))
    if (inherits(lok,"try-error")) return()
    X <- retrieve(destfile,lon=lon,lat=lat,verbose=verbose)
  }
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

## Function for generating urls to files available on a THREDDS data server (TDS). 
## Not sure if it works on other TDSs than the one here at metno.
## I had problems listing and downloading files from the same location on the metno TDS  
## and solved it by including two different base URLs, url.base and url.download. 
## Perhaps this can be done in an easier and more elegant way.
thredds.urls <- function(url.rel="raw/tas",pattern=".*EUR-11.*.nc",select=NULL,
                         url.base="http://thredds.met.no/thredds/catalog/postclim/data/CORDEX-EUR11",
                         url.download="http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11",
                         verbose=FALSE,...) {
  if(verbose) print("thredds.urls")
  url <- url.base
  if(!is.null(url.rel)) url <- paste(url,url.rel,sep="/")
  if(is.null(url.download)) url.download <- url.base
  if(!grepl(url,"catalog.html")) url <- paste(url,"catalog.html",sep="/")
  continue <- TRUE
  url.files <- NULL
  while(continue) {
    for(u in url) {
      txt <- readLines(u)
      txt <- txt[grep("href",txt)]
      files <- txt[grep(pattern,txt)]
      if(length(files)>0) {
        files <- gsub(paste("'>.*",sep=""),"",files)
        files <- gsub(paste(".*./",sep=""),"",files)
        u.data <- gsub(url.base,url.download,u)
        u.data <- gsub("catalog.html","",u.data)
        url.files <- c(url.files,paste(u.data,files,sep=""))
      }
      folders <- txt[grep("Folder",txt)]
      if(length(folders)>1) {
        folders <- gsub(".*href='|'>.*","",folders)
        url <- NULL
        for(i in seq(2,length(folders))) {
          u.i <- gsub("catalog.html",folders[i],u)
          url <- c(url,u.i)
        } 
      } else {
        continue <- FALSE
      }
    }
  }
  if(!is.null(select)) url.files <- url.files[select]
  return(url.files)
}

  #CORDEXList <- c("http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19500101-20051231.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_IPSL-IPSL-CM5A-MR_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v1_day_19510101-20051231.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_ICHEC-EC-EARTH_historical_r1i1p1_KNMI-RACMO22E_v1_day_19500101-20051231.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051230.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/tas/tas_EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19500101-20051231.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_IPSL-IPSL-CM5A-MR_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v1_day_19510101-20051231.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_day_19700101-20051231.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_ICHEC-EC-EARTH_historical_r1i1p1_KNMI-RACMO22E_v1_day_19500101-20051231.nc", "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_day_19700101-20051230.nc", 
  #                "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/raw/pr/pr_EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-20051231.nc")
  #
  #CORDEXAdjList <- c("http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/adjusted/tas/hist/tasAdjust_EUR-11_ICHEC-EC-EARTH_rcp45_r12i1p1_SMHI-RCA4_v1-METNO-QMAP-MESAN-1989-2010_day_19700101-20051231.nc", 
  #                   "http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11/adjusted/pr/hist/prAdjust_EUR-11_ICHEC-EC-EARTH_rcp45_r12i1p1_SMHI-RCA4_v1-METNO-QMAP-MESAN-1989-2010_day_19700101-20051231.nc")
