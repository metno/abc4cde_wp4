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
getCM <- function(url=NULL,destfile='CM.nc',lon=NULL,lat=NULL,force=FALSE,verbose=FALSE) {
  if(verbose) print("getCM")
  ## Retrieves the data
  if(is.null(url)) url <-
    'https://climexp.knmi.nl/CMIP5/monthly/tas/tas_Amon_ACCESS1-0_historical_000.nc'
  if (!file.exists(destfile)|force) lok <- try(download.file(url=url, destfile=destfile)) 
  if (inherits(loc,"try-error")) return()
  X <- retrieve(destfile,lon=lon,lat=lat)
  ## Collect information stored in the netCDF header
  cid <- getatt(destfile)
  ## Extract a time series for the area mean for 
  cid$area.mean <- aggregate.area(X,FUN='mean')
  cid$area.sd <- aggregate.area(X,FUN='sd')
  cid$url <- url
  cid$dates <- paste(range(index(X)),collapse=",")
  ncid <- nc_open(destfile)
  ## Collect information stored as model attributes
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
  if(is.null(destfile)) destfile <- paste(rep('GCM',length(select)),select,'.nc',sep='')
  ## Get the urls
  url <- cmip5.urls(varid=varid)[select] ## Get the URLs of the 
  ## Set up a list variable to contain all the metadata in sub-lists.
  X <- list()
  for (i in seq_along(select)) X[[paste('gcm',select[i],sep='.')]] <- getCM(url=url[i],destfile=destfile[i])
  return(X)
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
  if(is.null(destfile)) destfile <- paste(rep('CM',length(select)),select,'.nc',sep='')
  ## Get the urls
  url <- cordex.urls()[select]
  ## Set up a list variable to contain all the metadata
  X <- list()
  for (i in seq_along(select)) X[[paste('rcm',select[i],sep='.')]] <- getCM(url=url[i],destfile=destfile[i])
  return(X)
}

## Compute the common EOFs for GCMs and save the results for the front-end
commonEOFS.gcm <- function(select=1:9,varid='tas',destfile=NULL,
                           it='annual',is=NULL,verbose=FALSE) {
  if(verbose) print("commonEOFS.gcm")
  if(is.null(destfile)) destfile <- paste('GCM',select,'.nc',sep='')
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
  class(Z) <- c('dsensemble','eof','list')
  ceof <- Z
  save(ceof,file=paste('ceof.gcm',varid,it,'rda',sep='.'))
}

## Compute the common EOFs for RCMs save the results for the front-end
commonEOFS.rcm <- function(select=1:9,varid='tas',destfile=NULL,
                           it='annual',is=NULL,verbose=FALSE) {
  if(verbose) print("commonEOFS.rcm")
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
    xi <- attr(ceof,paste('appendix.',i,sep=''))
    rcmnames <- c(rcmnames,attr(xi,'model_id'))
    Z[[paste('rcm.',i+1,sep='')]] <- zoo(coredata(xi),order.by=index(xi))
    clim[[paste('rcm.',i+1,sep='')]] <- map.field(attr(X,paste('appendix.',i,sep="")))
  }
  attr(Z,'mean') <- clim
  attr(Z,'model_id') <- rcmnames
  class(Z) <- c('dsensemble','eof','list')
  ceof <- Z
  ##browser()
  save(ceof,file=paste('ceof.rcm',varid,it,'rda',sep='.'))
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
        
        ## Create output directory for the climate experiment
        #path.exp <- file.path(path,experiment[grep(iexp,experiment)],
        #                      fsep = .Platform$file.sep)
        #if (!file.exists(path.exp)) dir.create(path.exp)
        #if (verbose) print(path.exp[grep(iexp,experiment)])
        ## Define the output file
        #destfile <- paste(path.exp,varid[grep(ivar,varid)],sep="/") 
        #destfile <- paste(destfile,"_Amon_ens_",sep="")
        #destfile <- paste(destfile,iexp,sep="")
        #destfile <- paste(destfile,run.id,sep="_")
        #destfile <- paste(destfile,".nc",sep="")
        
        #if (!file.exists(destfile) | force) {
          ## Update output filename with attributes:
          urlfile  <- paste(url,ivar,sep="")             # add var directory
          urlfile  <- paste(urlfile,ivar,sep="/")        # add v.name
          urlfile  <- paste(urlfile,"_Amon_ens_",sep="") # add text
          urlfile  <- paste(urlfile,iexp,sep="")         # add exp.name
          urlfile  <- paste(urlfile,run.id,sep="_")      # add exp ID number
          urlfile  <- paste(urlfile,".nc",sep="")        # add file ext
        #}
        urlfiles <- c(urlfiles,urlfile)
        if (verbose) print(urlfile)
      }
      
    } # End for   
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
        
        ## Create output directory for the climate experiment
        #path.exp <- file.path(path,experiment[grep(iexp,experiment)],
        #                      fsep = .Platform$file.sep)
        #if (!file.exists(path.exp)) dir.create(path.exp)
        #if (verbose) print(path.exp[grep(iexp,experiment)])
        ## Define the output file
        #destfile <- paste(path.exp,varid[grep(ivar,varid)],sep="/") 
        #destfile <- paste(destfile,"EUR-44_cordex",sep="_")
        #destfile <- paste(destfile,iexp,"mon",sep="_")
        #destfile <- paste(destfile,run.id,sep="_")
        #destfile <- paste(destfile,".nc",sep="")
        
        #if (!file.exists(destfile) | force) {
          ## Update output filename with attributes:
          urlfile  <- paste(url,ivar,sep="")             # add var directory
          urlfile  <- paste(urlfile,ivar,sep="/")        # add v.name
          urlfile  <- paste(urlfile,"EUR-44_cordex",sep="_") # add text
          urlfile  <- paste(urlfile,iexp,"mon",sep="_")         # add exp.name
          urlfile  <- paste(urlfile,run.id,sep="_")      # add exp ID number
          urlfile  <- paste(urlfile,".nc",sep="")        # add file ext
        #}
        urlfiles <- c(urlfiles,urlfile)
        if (verbose) print(urlfile)
      }
      
    } # End for   
  }
  return(urlfiles[-1])
}

metaextract <- function(x=NULL, verbose=FALSE) {
  if(verbose) print("metaextract")
  if (is.null(x)) x <- c(getGCMs(verbose=verbose),getRCMs(verbose=verbose))
  gcms <- names(x)
  n <- length(gcms)
  if(verbose) print(gcms)
  for(i in seq_along(x)) {
    xx <- x[[gcms[i]]]
    if(is.null(xx$project_id)) {
      print(paste("Warning! project_id is not specified in",xx$filename))
      yi <- NULL
    } else if(grepl("cmip",tolower(xx$project_id))) {
      yi <- metaextract.cmip(xx,verbose=verbose)
    } else if(grepl("cordex",tolower(xx$project_id))) {
      yi <- metaextract.cordex(xx,verbose=verbose)
    }
    print(i)
    if(i==1) {
      Y <- matrix(NA,ncol=ncol(yi),nrow=n)
      colnames(Y) <- colnames(yi)
      Y[i,] <- yi
    } else {
      cn.all <- unique(c(colnames(Y),colnames(yi)))
      Y.new <- matrix(NA,ncol=length(cn.all),nrow=n)
      colnames(Y.new) <- cn.all
      j <- sapply(colnames(Y),function(x) which(cn.all==x))
      Y.new[1:(i-1),j] <- Y[1:(i-1),]
      #for(k in 1:(i-1)) Y.new[k,j] <- Y[k,]
      for(cn in colnames(yi)) {
        Y.new[i,colnames(Y.new)==cn] <- yi[colnames(yi)==cn]
      }
      Y <- Y.new
    }
  }
  Y -> metaextract
  save(metaextract,file='metaextract.rda')
  return(Y)
}

metaextract.cmip <- function(x=NULL, verbose=FALSE) {
  if(verbose) print("metaextract.cmip")
  ## argument 'x' is input from getGCMs, getRCMs, testGCM, etc
  if (is.null(x)) x <- c(getGCMs(verbose=verbose),getRCMs(verbose=verbose))
  
  if(!inherits(x,"list")) x <- list(gcm.1=x)
  gcms <- names(x)
  n <- length(gcms)
  
  for (i in 1:n) {
    xx <- x[[gcms[i]]]
    project_id <- NA; url <- NA; filename <- NA; dim <- NA; dates <- NA
    var <- NA; longname <- NA; vunit <- NA; vid <- NA
    res <- NA; lon.rng <- NA; lon.unit <- NA; lat.rng <- NA; lat.unit <- NA
    experiment_id <- NA; frequency <- NA; creation_date <- NA; tracking_id <- NA
    gcm <- NA; gcm.rip <- NA; gcm.v <- NA; gcm.realm <- NA
    if(!is.null(xx$dim)) dim <- paste(names(xx$dim),collapse=",")
    if(!is.null(names(xx$var))) {
      var <- names(xx$var)#[!grepl("time",names(xx$var))]
      if(!is.null(xx$var[[1]]$longname)) longname <- sapply(var, function(x) xx$var[[x]]$longname)
      if(!is.null(xx$var[[1]]$units)) vunit <- sapply(var, function(x) xx$var[[x]]$units)
      if(!is.null(xx$var[[1]]$id$id)) vid <- sapply(var, function(x) xx$var[[x]]$id$id)
    }
    if(!is.null(names(xx$dim))) {
      if(!is.null(xx$dim$lat$vals)) {
        res <- diff(xx$dim$lat$vals)[1]
        lat.rng <- paste(range(xx$dim$lat$vals),collapse=",")
      }
      if(!is.null(xx$dim$lon$vals)) lon.rng <- paste(range(xx$dim$lon$vals),collapse=",")
      if(!is.null(xx$dim$lat$units)) lat.unit <- xx$dim$lat$units
      if(!is.null(xx$dim$lon$units)) lon.unit <- xx$dim$lon$units
    }
    for(mi in c("url","filename","dates")) {
      if(!is.null(xx[[mi]])) eval(parse(text=paste(mi," <- xx$",mi,sep="")))
    }
    for(mi in c("project_id","experiment_id","frequency","creation_date","tracking_id")) {
      if(!is.null(xx$model[[mi]])) eval(parse(text=paste(mi," <- xx$model$",mi,sep="")))
    }
    if(!is.null(xx$model$model_id)) gcm <- xx$model$model_id
    if(!is.null(xx$model$parent_experiment_rip)) gcm.rip <- xx$model$parent_experiment_rip
    if(!is.null(xx$model$version_number)) gcm.v <- xx$model$version_number
    if(!is.null(xx$model$modeling_realm)) gcm.realm <- xx$model$modeling_realm
    mx <- data.frame(project_id=project_id, url=url, filename=filename,
                     dim=paste(dim,collapse=","), dates=dates, var=paste(var,collapse=","),
                     longname=paste(longname,collapse=","), unit=paste(vunit,collapse=","),
                     var_id=paste(vid,collapse=","), 
                     resolution=res, lon=lon.rng, lon_unit=lon.unit, lat=lat.rng, lat_unit=lat.unit,
                     experiment_id=experiment_id, frequency=frequency, 
                     creation_date=creation_date, tracking_id=tracking_id,
                     gcm=gcm, gcm_rip=gcm.rip, gcm_version=gcm.v, gcm_realm=gcm.realm)
    meta <- names(mx)
    m <- length(meta)
    if (i==1) {
      X <- matrix(rep("NA",n*m),n,m) ## set up a matrix
      colnames(X) <- meta; rownames(X) <- gcms
    }
    for (ii in 1:m) {
      if(!is.na(mx[[meta[ii]]])) {
        y <- as.character(mx[[meta[ii]]])
        X[i,ii] <- y
      }
    }
  }
  return(X)
}

metaextract.cordex <- function(x=NULL, verbose=FALSE) {
  if(verbose) print("metaextract.cordex")
  ## argument 'x' is input from getGCMs, getRCMs, testGCM, etc
  if (is.null(x)) x <- c(getGCMs(),getRCMs())
  
  if(!inherits(x,"list")) x <- list(gcm.1=x)
  gcms <- names(x)
  n <- length(gcms)
  
  for (i in 1:n) {
    xx <- x[[gcms[i]]]
    project_id <- NA; url <- NA; filename <- NA; dim <- NA; dates <- NA
    var <- NA; longname <- NA; vunit <- NA; vid <- NA
    res <- NA; lon.rng <- NA; lon.unit <- NA; lat.rng <- NA; lat.unit <- NA
    experiment_id <- NA; frequency <- NA; creation_date <- NA; tracking_id <- NA
    gcm <- NA; gcm.rip <- NA; gcm.v <- NA; gcm.realm <- NA
    if(!is.null(xx$dim)) dim <- paste(names(xx$dim),collapse=",")
    if(!is.null(names(xx$var))) {
      var <- names(xx$var)#[!grepl("time",names(xx$var))]
      if(!is.null(xx$var[[1]]$longname)) longname <- sapply(var, function(x) xx$var[[x]]$longname)
      if(!is.null(xx$var[[1]]$units)) vunit <- sapply(var, function(x) xx$var[[x]]$units)
      if(!is.null(xx$var[[1]]$id$id)) vid <- sapply(var, function(x) xx$var[[x]]$id$id)
    }
    if(!is.null(names(xx$dim))) {
      if(!is.null(xx$dim$lat$vals)) {
        res <- diff(xx$dim$lat$vals)[1]
        lat.rng <- paste(range(xx$dim$lat$vals),collapse=",")
      }
      if(!is.null(xx$dim$lon$vals)) lon.rng <- paste(range(xx$dim$lon$vals),collapse=",")
      if(!is.null(xx$dim$lat$units)) lat.unit <- xx$dim$lat$units
      if(!is.null(xx$dim$lon$units)) lon.unit <- xx$dim$lon$units
    }
    for(mi in c("url","filename","dates")) {
      if(!is.null(xx[[mi]])) eval(parse(text=paste(mi," <- xx$",mi,sep="")))
    }
    for(mi in c("project_id","experiment_id","frequency","creation_date","tracking_id")) {
      if(!is.null(xx$model[[mi]])) eval(parse(text=paste(mi," <- xx$model$",mi,sep="")))
    }
    if(!is.null(xx$model$driving_model_id)) gcm <- xx$model$driving_model_id
    if(!is.null(xx$model$driving_model_ensemble_member)) gcm.rip <-
                                         xx$model$driving_model_ensemble_member
    if(!is.null(xx$model$model_id)) rcm <- xx$model$model_id
    if(!is.null(xx$model$CORDEX_domain)) rcm.domain <- xx$model$CORDEX_domain
    if(!is.null(xx$model$rcm_version_id)) rcm.v <- xx$model$rcm_version_id
    mx <- data.frame(project_id=project_id, url=url, filename=filename,
                     dim=paste(dim,collapse=","), dates=dates, var=paste(var,collapse=","),
                     longname=paste(longname,collapse=","), unit=paste(vunit,collapse=","),
                     var_id=paste(vid,collapse=","), 
                     resolution=res, lon=lon.rng, lon_unit=lon.unit, lat=lat.rng, lat_unit=lat.unit,
                     experiment_id=experiment_id, frequency=frequency, 
                     creation_date=creation_date, tracking_id=tracking_id,
                     gcm=gcm, gcm_rip=gcm.rip, rcm=rcm, rcm_domain=rcm.domain, rcm_version=rcm.v)
    meta <- names(mx)
    m <- length(meta)
    if (i==1) {
      X <- matrix(rep("NA",n*m),n,m) ## set up a matrix
      colnames(X) <- meta; rownames(X) <- gcms
    }
    for (ii in 1:m) {
      if(!is.na(mx[[meta[ii]]])) {
        y <- as.character(mx[[meta[ii]]])
        X[i,ii] <- y
      }
    }
  }
  return(X)
}

