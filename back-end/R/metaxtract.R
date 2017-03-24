## ABC4CDE/DECM WP4 R-scripts for propotype climate data storage (CDS) tool for 
## processing data, metadata and evaluationand quality control (EQC)
## Rasmus.Benestad@met.no 2017-02-15

## Function that extracts metadata stored in lists and sublists so that it's easier to search and
## organise into tables/spreadsheets. 

metaxtract <- function(x=NULL) {
  ## argument 'x' is input from getGCMs, getRCMs, testGCM, etc
  if (is.null(x)) x <- c(getGCMs(),getRCMs())
  
  gcms <- names(x)
  n <- length(gcms)
  print(gcms)
  
  for (i in 1:n) {
    ## (the id (xx$id) tag is not a unique identifier of the specific model or file)
    xx <- x[[gcms[i]]] ## Tidier - use xx for the next level in the list structure
    project_id <- NA; url <- NA; filename <- NA; dim <- NA; dates <- NA
    var <- NA; longname <- NA; vunit <- NA; vid <- NA
    res <- NA; lon.rng <- NA; lon.unit <- NA; lat.rng <- NA; lat.unit <- NA
    experiment_id <- NA; frequency <- NA; creation_date <- NA; tracking_id <- NA
    gcm <- NA; gcm.rip <- NA; gcm.v <- NA; gcm.realm <- NA
    rcm <- NA; rcm.domain <- NA; rcm.v <- NA
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
    if(grepl("cordex",tolower(project_id))) {
      if(!is.null(xx$model$driving_model_id)) gcm <- xx$model$driving_model_id
      if(!is.null(xx$model$driving_model_ensemble_member)) gcm.rip <- xx$model$driving_model_ensemble_member
      if(!is.null(xx$model$model_id)) rcm <- xx$model$model_id
      if(!is.null(xx$model$CORDEX_domain)) rcm.domain <- xx$model$CORDEX_domain
      if(!is.null(xx$model$rcm_version_id)) rcm.v <- xx$model$rcm_version_id
    } else if(grepl("cmip",tolower(project_id))) {
      if(!is.null(xx$model$model_id)) gcm <- xx$model$model_id
      if(!is.null(xx$model$parent_experiment_rip)) gcm.rip <- xx$model$parent_experiment_rip
      if(!is.null(xx$model$version_number)) gcm.v <- xx$model$version_number
      if(!is.null(xx$model$modeling_realm)) gcm.realm <- xx$model$modeling_realm
    } else {
      print("Unexpected project id:",mx$project_id)
    }
    mx <- data.frame(project_id=project_id, url=url, filename=filename,
                     dim=paste(dim,collapse=","), dates=dates, var=paste(var,collapse=","),
                     longname=paste(longname,collapse=","), unit=paste(vunit,collapse=","),
                     var_id=paste(vid,collapse=","), 
                     resolution=res, lon=lon.rng, lon_unit=lon.unit, lat=lat.rng, lat_unit=lat.unit,
                     experiment_id=experiment_id, frequency=frequency, 
                     creation_date=creation_date, tracking_id=tracking_id,
                     gcm=gcm, gcm_rip=gcm.rip, gcm_version=gcm.v, gcm_realm=gcm.realm,
                     rcm=rcm, rcm_domain=rcm.domain, rcm_version=rcm.v)
    meta <- names(mx)
    m <- length(meta)
    if (i==1) {
      X <- matrix(rep("NA",n*m),n,m) ## set up a matrix
      colnames(X) <- meta; rownames(X) <- gcms
    }
    for (ii in 1:m) {
      if(!is.na(mx[[meta[ii]]])) {
        y <- as.character(mx[[meta[ii]]])
        #print(paste(i,ii,meta[[ii]],y))
        X[i,ii] <- y
      }
    }
  }
  X -> metaxtract
  save(metaxtract,file='metaxtract.rda')
  return(X)
}