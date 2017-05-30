## ABC4CDE/DEMC - R-script for prototype tool WP4
## andreas.dobler@met.no  Oslo, Norway, 2017-05-16
##

##  Function of an R-package that retrieves metadata from files located on OPeNDAP servers.
require(PCICt) ##Used to get the time axes right for 360 day and non-leapyear calendars.

#Function to extract the metadata form NetCDF files on OPeNDAP servers
metaextract.opendap <- function(url=NULL, verbose=FALSE) {
  if(verbose) print("metaextract_opendap")
  
  #Make connection to file and get global attributes, variable and dimension lists
  nc <- nc_open(url)
  globat <- ncatt_get(nc,0)
  varlist <- nc$var
  dimlist <- nc$dim
  nc_close(nc)
  
  #Get the time range in the data (and convert to standard calender)
  tunit <- dimlist$time$units
  tsplit <- unlist(strsplit(tunit,split=" "))
  torigin <- paste(tsplit[3:length(tsplit)],collapse=" ")
  tvals <- dimlist$time$vals
  tdiff <- as.difftime(tvals,units=tsplit[1])
  tdiffsecs <- as.numeric(tdiff, units = "secs")
  caltype <- dimlist$time$calendar
  timeline <- as.PCICt(torigin,cal=caltype) + tdiffsecs
  
  #Return a list with: Global attributes, list of variables, list of dimensions, time range
  return(list(globat=globat,varlist=varlist,dimlist=dimlist,url=url,range=range(timeline)))
}

## Function to extract the metadata from local NetCDF files
## [Change so that the opendap alternative can be used for data from the KNMI explorer too]
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
      yi <- metaextract.cmip5(xx,verbose=verbose)
    } else if(grepl("cordex",tolower(xx$project_id))) {
      yi <- metaextract.cordex(xx,verbose=verbose)
    }
    if(verbose) print(i)
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
      for(cn in colnames(yi)) {
        Y.new[i,colnames(Y.new)==cn] <- yi[colnames(yi)==cn]
      }
      Y <- Y.new
    }
  }
  Y -> meta
  save(meta,file='metaextract.rda')
  return(Y)
}

metaextract.cmip5 <- function(x=NULL, verbose=FALSE) {
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
    for(mi in c("url","filename","dates","frequency")) {
      if(!is.null(xx[[mi]])) eval(parse(text=paste(mi," <- xx$",mi,sep="")))
    }
    for(mi in c("project_id","experiment_id","creation_date","tracking_id")) {
      if(!is.null(xx$model[[mi]])) eval(parse(text=paste(mi," <- xx$model$",mi,sep="")))
    }
    if(!is.null(xx$model$model_id)) gcm <- xx$model$model_id
    if(!is.null(xx$model$parent_experiment_rip)) gcm.rip <- xx$model$parent_experiment_rip
    if(!is.null(xx$model$version_number)) gcm.v <- xx$model$version_number
    if(!is.null(xx$model$modeling_realm)) gcm.realm <- xx$model$modeling_realm
    mx <- data.frame(project_id=project_id, url=url, filename=filename,
                     dim=paste(dim,collapse=","), dates=dates, var=paste(var,collapse=","),
                     longname=paste(longname,collapse=","), unit=paste(vunit,collapse=","),
                     #var_id=paste(vid,collapse=","), 
                     resolution=res, lon=lon.rng, lon_unit=lon.unit, lat=lat.rng, lat_unit=lat.unit,
                     experiment_id=experiment_id, frequency=frequency, 
                     creation_date=creation_date, #tracking_id=tracking_id,
                     gcm=gcm, gcm_rip=gcm.rip)#, gcm_version=gcm.v, gcm_realm=gcm.realm)
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
    if(!is.null(xx$dim)) dim <- paste(names(xx$dim),collapse=",")[!grepl("bnds",names(xx$var))]
    if(!is.null(names(xx$var))) {
      var <- names(xx$var)[!grepl("bnds",names(xx$var))]
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
                     #var_id=paste(vid,collapse=","), 
                     resolution=res, lon=lon.rng, lon_unit=lon.unit, lat=lat.rng, lat_unit=lat.unit,
                     experiment_id=experiment_id, frequency=frequency, 
                     creation_date=creation_date, #tracking_id=tracking_id,
                     gcm=gcm, gcm_rip=gcm.rip, rcm=rcm)#, rcm_domain=rcm.domain, rcm_version=rcm.v)
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


