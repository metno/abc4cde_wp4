## ABC4CDE/DEMC - R-script for prototype tool WP4
## Rasmus.Benestad@met.no  Oslo, Norway, 2017-02-14
##

require(esd) ## This code builds on the esd package: https://github.com/metno/esd
require(raster) 
require(rgdal)

##  Function of an R-package that retrieves data files with CMIP or CORDEX results

#Helper function to find files linux environment (add Windows counterpart)
find.file <- function(filename){
  command <- paste("find $HOME -name",filename,sep=" ")
  fullpath <- system(command,intern=T)
  return(fullpath)
}

#apply mask to a zoo object by setting values outside the mask to NA
mask.zoo <- function(zoo.object,mask){
  mask <- flip(mask,direction='y')
  zoo.object[,which(is.na(getValues(mask)))] <- NA
  return(zoo.object)
}

#Search and read a shapefile
get.shapefile <- function(filename=NULL,with.path=F){
  fullname <- filename
  if(!with.path){
    fullname <- find.file(filename)
  }
  readOGR(fullname,verbose=F)
}

#Apply a set of cdo commands on a grib/netcdf file.
#Several commands can be piped.
cdo.command <- function(commands,input,infile,outfile){
  cdo.coms <- array()
  separators <- array(" ",dim=length(commands))
  separators[which(is.na(match(input,"")))] <- ","
  for(i in 1:length(separators)){
    cdo.coms[i]  <- paste(commands[i],input[i],sep=separators[i])
  }
  system.command <- paste("cdo",paste(cdo.coms,collapse=" "),infile,outfile,sep=" ")
  system(system.command,wait=T)
}

#Unzip a gz package
gunzip <- function(filename){
  system.command <- paste("gunzip",filename)
  system(system.command,wait=T)
}


#Get grid boxes belonging to a SREX region and calculate some basic statistics for it.
get.srex.region <- function(destfile,region=NULL,print.srex=F,verbose=F){ 
  home <- system("echo $HOME",intern=T)
  shape <-  get.shapefile("referenceRegions.shp")
  X <- retrieve(destfile,lon=NULL,lat=NULL,verbose=verbose)
  srex <- list()
  if(is.null(region)){
    for (i in 1:length(levels(shape$LAB))){
      polygon <- shape[i,]
      mask <- gen.mask.srex(destfile=destfile, mask=polygon, ind=F, inverse=F, mask.values=1)
      if(verbose){
        if(i==1){
          plot(shape)
        }
        plot.mask <- mask
        extent(plot.mask) <- c(-180,180,-90,90)
        projection(plot.mask) <- projection(shape)
        plot(plot.mask,col=rainbow(100, alpha=0.35)[sample(1:100,1)],legend=F,add=T)
      }
      name <- levels(shape$NAME)[i]
      X.region <- mask.zoo(X,mask)
      srex[[name]]$name <- name
      srex[[name]]$label <- levels(shape$LAB)[i]
      srex[[name]]$area.mean <- aggregate.area(X.region,FUN="mean",na.rm=T)
      srex[[name]]$area.sd <- aggregate.area(X.region,FUN="sd",na.rm=T)
    }  
  }else{
    polygon <- shape[levels(shape$LAB)==region,]
    mask <- gen.mask.srex(destfile=destfile, mask=polygon, ind=F, inverse=F, mask.values=1)
   if(verbose){
      plot(shape)
      plot.mask <- mask
      extent(plot.mask) <- c(-180,180,-90,90)
      projection(plot.mask) <- projection(shape)
      plot(plot.mask,col=rainbow(100, alpha=0.35)[sample(1:100,1)],legend=F,add=T)
   }
    name <- levels(shape$NAME)[i]
    X.region <- mask.zoo(X,mask)
    srex[[name]]$name <- name
    srex[[name]]$label <- levels(shape$LAB)[i]
    srex[[name]]$area.mean <- aggregate.area(X.region,FUN="mean",na.rm=T)
    srex[[name]]$area.sd <- aggregate.area(X.region,FUN="sd",na.rm=T) 
  }

  if(print.srex){
    print("Region names in alphabetical order and the corresponding label to be used when selecting the region:")
    print(data.frame(NAME=gsub("\\[[^\\]]*\\]", "", levels(shape$NAME), perl=TRUE),
                        LABEL=levels(shape$LAB)))
    return()
  }
  return(srex)
}

#Create a raster mask for the selected SREX sub-region from the CMIP5 netcdf file.
gen.mask.srex <- function(destfile, mask=NULL, ind=F, inverse=F, mask.values=1){
  print(destfile)
  r <- raster(destfile)
  r <- setValues(r,NA)
  extent.r <- extent(r)
  if(extent.r[2]==360) extent(r) <- c(-180,180,-90,90)
  indices <- extract(r,mask,cellnumbers=T)[[1]][,1]
  if(extent(mask)[2]>180){
    extent(r) <- c(180,540,-90,90)
  }
  indices <- sort(c(indices,extract(r,mask,cellnumbers=T)[[1]][,1]))
  if(inverse){
    tmp <- seq(1,length(getValues(r)))
    indices <- tmp[which(is.na(match(tmp,indices)))]
  }
  mask.raster <- r
  extent(mask.raster) <- c(0,360,-90,90)
  mask.raster[indices] <- mask.values
  if(ind)return(indices)
  return(mask.raster)
}

getatt <- function(fname) {
  ## Reads and extracts the attribute information in a netCDF files and stores this in a list object## 
  ## This is part of the information stored in the metadatabase
  ncid <- nc_open(fname)
  nc_close(ncid)
  return(ncid)
}

#Call to a python script which downloads data from the public ECMWF data server
python.getEra <- function(start,end,variable,steps,type,stream,outfile){
  script <- "python ~/abc4cde_wp4/back-end/python/getMonthlyERA.py"
  system.command <- paste(script," -f ",start," -l ",end," -v ",variable,
                          " -s ",steps," -t ",type," -r ",stream," -o ",outfile, sep="")
  system(system.command,wait=T)
}

#Retrieve monthly data for 2m temperature and precipitation from the ECMWF public repository.
#The use of this function requires that the ECMWF key and python libraries are installed on the machine.
#See instructions in https://software.ecmwf.int/wiki/display/WEBAPI/Access+ECMWF+Public+Datasets
#The function also requires that cdo is installed on the operating computer.
getERA <- function(variable.name,start=1979,end=2016,griddes="cmip_1.25deg_to_2.5deg.txt",destfile=NULL){
  griddes <- find.file(griddes)
  if(any(match(c("tas","tmp","temp","temperature","t2m"),variable.name,nomatch=0))){
    varID<- "167.128"
    stream <- "moda"
    steps <- "0"
    type <- "an"
    commands <- c("-f","nc","-copy","-remapcon","-chname")
    input <- c("","","",griddes,"2t,tas")
#    seps <- c(" ","",",",",")
  }else if(any(match(c("pre","prc","prec","precipitation","pr"),variable.name,nomatch=0))){
    varID<- "228.128"
    stream <- "mdfa"
    steps <- "12-24/24-36" # Select the 24 and 36h 12-hour long forecast in order to reduce the spin-up effect.
    type <- "fc"
    commands <- c("-f","nc","-copy","-monsum","-remapcon","-chname")
    input <- c("","","","",griddes,"2t,tas")
#    seps <- c(" ","","",",",",")
  }
  griddes <- find.file(griddes)
  if(is.null(destfile)) destfile <- paste("era-interim_monthly_",paste(start,end,sep="-"),"_",variable.name,".grib",sep="")
  if(!file.exists(destfile))python.getEra(start, end, varID, steps, type, stream, destfile)
  outfile <- paste(gsub('.{5}$', '',destfile),"2.5deg",'nc',sep=".")
  if(!file.exists(outfile)) cdo.command(commands,input,infile=destfile,outfile=outfile)
  X <- retrieve(outfile)
  cid <- getatt(outfile)
  cid$area.mean <- aggregate.area(X,FUN='mean')
  cid$area.sd <- aggregate.area(X,FUN='sd')
  cid$url <- NA
  cid$dates <- paste(range(index(X)),collapse=",")
  ncid <- nc_open(outfile)
  model <- ncatt_get(ncid,0)
  nc_close(ncid)
  cid$model <- model
  cid$project_id <- cid$model$project_id
  cid$srex <- get.srex.region(outfile,region=NULL,print.srex=F,verbose=F)
  return(cid)
}

#A test function to retrieve CRU data from CEDA databases.
getCRU <- function(username,passwd,variable="tmp",version="4.00",griddes="cmip_1.25deg_to_2.5deg.txt",destfile=NULL,time.resol=NULL){
  if(any(match(c("tas","tmp","temp","temperature","t2m"),variable))){
    variable <- "tmp"
  }else if(any(match(c("pre","prc","prec","precipitation","pr"),variable))){
    variable <- "pre"
  }
  cert <- paste(username,passwd,sep=":")
  url <- paste("ftp.ceda.ac.uk/badc/cru/data/cru_ts",paste("cru_ts",version,sep="_"),"data",variable,sep="/")
  if(is.null(destfile)) destfile <- paste(paste("cru_ts",version,sep=""),"1901.2015",variable,"dat.nc.gz",sep=".")
  if(!file.exists(destfile))try(download.file(url=paste(paste("ftp://",paste(cert,url,sep="@"),sep=""),destfile,sep="/"),
                                              destfile=destfile, mode="wb"))
  gunzip(destfile)
  destfile <- paste(gsub('.{5}$', '',destfile),"nc",sep="")
  outfile <- paste(gsub('.{5}$', '',destfile),"2.5deg.",'nc',sep="")
  if(!file.exists(outfile)){
    griddes <- find.file(griddes)
    commands <- c("-f","-copy","-remapcon")
    input <- c("nc","",griddes)
    seps <- c(" ","",",")
    cdo.command(commands,input,seps,infile=destfile,outfile=outfile) 
  }
  X <- retrieve(outfile)
  cid <- getatt(outfile)
  cid$area.mean <- aggregate.area(X,FUN='mean')
  cid$area.sd <- aggregate.area(X,FUN='sd')
  cid$url <- NA
  cid$dates <- paste(range(index(X)),collapse=",")
  ncid <- nc_open(outfile)
  model <- ncatt_get(ncid,0)
  nc_close(ncid)
  cid$model <- model
  cid$project_id <- cid$model$project_id
  return(cid)
}

#Get monthly CFSR data and interpolate it to common 2.5 degree grid.
getCFSR <- function(variable="t2m",destfile=NULL,lon=NULL,lat=NULL,verbose=T,griddes="cmip_1.25deg_to_2.5deg.txt"){
  url.path <- "https://climexp.knmi.nl/CFSR"
  griddes <- find.file(griddes)
  if(variable=="tas"){
    filename <- "cfsr_tmp2m.nc"
    commands <- c("-f","-copy","-remapcon","-monavg","-chname")
    input <- c("nc","",griddes,"","TMP_2maboveground,t2m")
    seps <- c(" ","",","," ",",")
  }else if(variable=="pr"){
    filename <- "cfsr_prate.nc"
    commands <- c("-f","-copy","-remapcon","-monavg","-chname")
    input <- c("nc","",griddes,"","PRATE_surface,Pr")
    seps <- c(" ","",",","",",")
  }
  if(!file.exists(filename))download.file(paste(url.path,filename,sep="/"),destfile=filename)
  if(is.null(destfile)) destfile <- paste(sub("\\.[[:alnum:]]+$", "", filename, perl=TRUE),"mon.nc",sep="_")
  if(!file.exists(destfile)) cdo.command(commands,input,seps,pipe=T,infile=filename,outfile=destfile)
  X <- retrieve(destfile,lon=lon,lat=lat,verbose=verbose)
  cid <- getatt(destfile) 
  cid$url <- paste(url.path,filename,sep="/")
  cid$area.mean <- aggregate.area(X,FUN='mean',na.rm=T)
  cid$area.sd <- aggregate.area(X,FUN='sd',na.rm=T)
  ncid <- nc_open(destfile)
  model <- ncatt_get(ncid,0)
  nc_close(ncid)
  cid$model <- model
  cid$srex <- get.srex.region(destfile,region=NULL,print.srex=F,verbose=F)
  #  file.remove(filename)
  return(cid)
}

#Get daily EOBS data and convert it to monthly averages
getEOBS <- function(variable="tas", destfile=NULL, resolution="0.50", version="14"){
  url.path <- "http://www.ecad.eu/download/ensembles/data/Grid_0.50deg_reg"
  if(variable=="tas"){
    filename <- "tg_0.50deg_reg_v14.0.nc.gz"
  }else if(variable=="pr"){
    filename <- "rr_0.50deg_reg_v14.0.nc.gz"
  }else{
    return("Not implemented yet!")
  }
  if(!file.exist(filename)) download.file(paste(url.path,filename,sep="/"),destfile=filename)
  gunzip(filename)
  filename <- sub("\\.[[:alnum:]]+$", "", filename, perl=TRUE)
  if(is.null(destfile)) destfile <- paste(sub("\\.[[:alnum:]]+$", "", filename, perl=TRUE),"mon.nc",sep="_")
  commands <- c("-f","-copy","-monavg")
  input <- c("nc","","")
  seps <- c(" ","","")
  if(!file.exist(destfile)) cdo.command(commands,input,seps,pipe=F,infile=filename,outfile=destfile)
  #  file.remove(filename)
  X <- retrieve(destfile,lon=lon,lat=lat,verbose=verbose)
  cid <- getatt(destfile) 
  cid$url <- paste(url.path,filename,sep="/")
  cid$area.mean <- aggregate.area(X,FUN='mean',na.rm=T)
  cid$area.sd <- aggregate.area(X,FUN='sd',na.rm=T)
  ncid <- nc_open(destfile)
  model <- ncatt_get(ncid,0)
  nc_close(ncid)
  cid$model <- model
  #  file.remove(filename)
  return(cid)
}

## Generic function to retrieve climate model (CM) file from the KNMI ClimateExplorer
getCM <- function(url=NULL,destfile='CM.nc',lon=NULL,lat=NULL,force=FALSE,verbose=FALSE) {
  if(verbose) print("getCM")
  ## Retrieves the data
  if(is.null(url)) url <-
    'https://climexp.knmi.nl/CMIP5/monthly/tas/tas_Amon_ACCESS1-0_historical_000.nc'
  if (!file.exists(destfile)|force){
    loc <- try(download.file(url=url, destfile=destfile))
    if (inherits(loc,"try-error")) return()
  } 
  X <- retrieve(destfile,lon=lon,lat=lat,verbose=verbose)
  ## Collect information stored in the netCDF header
  cid <- getatt(destfile)
  ## Extract a time series for the area mean 
  cid$area.mean <- aggregate.area(X,FUN='mean',na.rm=T)
  cid$area.sd <- aggregate.area(X,FUN='sd',na.rm=T)
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

getGCMs <- function(select=1:9,varid='tas',destfile=NULL,verbose=FALSE,get.srex=TRUE,region=NULL,print.srex=FALSE){
  if(verbose) print("getGCMs SREX regions")
  ## Set destfile
  if(is.null(destfile)) destfile <- paste(rep('GCM',length(select)),select,'.',varid,'.nc',sep='')
  url <- cmip5.urls(varid=varid)[select] ## Get the URLs of the 
  ## Set up a list variable to contain all the metadata in sub-lists.
  X <- list()
  for (i in seq_along(select)) {
    if(verbose) print(paste("Get gcm.",select[i],sep=''))
    X[[paste('gcm',varid,select[i],sep='.')]] <-
        getCM(url=url[i],destfile=destfile[i],verbose=verbose)
    if(get.srex)X[[paste('gcm',varid,select[i],sep='.')]]$srex <-
        get.srex.region(destfile=destfile[i],region=region,print.srex=print.srex,verbose=verbose)
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
getRCMs <- function(select=1:9,varid='tas',destfile=NULL,verbose=FALSE,region=NULL) {
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

## Compute the common EOFs for GCMs and save the results for the front-end
commonEOFS.gcm <- function(select=1:9,varid='tas',destfile=NULL,
                           it='annual',is=NULL,verbose=FALSE) {
  if(verbose) print("commonEOFS.gcm")
  if(is.null(destfile)) destfile <- paste('GCM',select,'.',varid,'.nc',sep='')
  getGCMs(select=select,varid=varid,destfile=destfile)
  X <- NULL
  for (fname in destfile) {
    if(verbose) print(paste("retrieve",fname))
    x <- retrieve(fname,verbose=verbose)
    if (!is.null(it)) {
      if (tolower(it)=='annual') x <- annual(x,verbose=verbose) else
                                 x <- subset(x,it=it,is=is,verbose=verbose)
    }
    if (is.null(X)) X <- x else X <- combine(X,x,verbose=verbose)
  }
  if(verbose) print("Calculate common EOF")
  ceof <- EOF(X,verbose=verbose)
  plot(ceof)
  
  ## Need to reformat the ceof-object to fit the set-up for the R-shiny app in the front-end.
  if(verbose) print("Reformat the common EOF object")
  x1 <- coredata(ceof); attributes(x1) <- NULL; dim(x1) <- dim(ceof)
  eof <- as.eof(ceof); attr(eof,"standard.error") <- NULL
  Z <- list(info='CMIP5 runs',eof=eof,rcm.1=zoo(x1,order.by=index(ceof)))
  clim <- list(rcm.1=map.field(X,plot=FALSE))
  gcmnames <- attr(ceof,'model_id')
  gcmrip <- attr(ceof,'parent_experiment_rip')
  for (i in 1:attr(ceof,'n.apps')) {
    x1 <- attr(ceof,paste('appendix.',i,sep=''))
    Z[[paste('rcm.',i+1,sep='')]] <- zoo(coredata(x1),order.by=index(x1))
    gcmnames <- c(gcmnames,attr(x1,'model_id'))
    gcmrip <- c(gcmrip,attr(x1,'parent_experiment_rip'))
    clim[[paste('rcm.',i+1,sep='')]] <- map.field(attr(X,paste('appendix.',i,sep='')))
  }
  attr(Z,'mean') <- clim
  attr(Z,'model_id') <- list(gcm=gcmnames,gcm_rip=gcmrip)
  class(Z) <- c('dsensemble','eof','list')
  ceof <- Z
  save(ceof,file=paste('ceof.gcm',varid,it,'rda',sep='.'))
  invisible(ceof)
}

## Compute the common EOFs for RCMs save the results for the front-end
commonEOFS.rcm <- function(select=1:9,varid='tas',destfile=NULL,
                           it='annual',is=NULL,verbose=FALSE) {
  if(verbose) print("commonEOFS.rcm")
  if(is.null(destfile)) destfile <- paste(rep('CM',length(select)),select,'.',varid,'.nc',sep='')
  getRCMs(select=select,varid=varid,destfile=destfile,verbose=verbose)
  X <- NULL
  for (fname in destfile) {
    if(verbose) print(paste("retrieve",fname))
    x <- retrieve(fname)
    if (!is.null(it)) {
      if (tolower(it)=='annual') x <- annual(x,verbose=verbose) else
        x <- subset(x,it=it,is=is,verbose=verbose)
    }
    if (is.null(X)) X <- x else X <- combine(X,x,verbose=verbose)
  }
  if(verbose) print("Calculate common EOF")
  ceof <- EOF(X,verbose=verbose)
  plot(ceof)
  
  ## Need to reformat the ceof-object to fit the set-up for the R-shiny app in the front-end.
  if(verbose) print("Reformat the common EOF object")
  x1 <- coredata(ceof); attributes(x1) <- NULL; dim(x1) <- dim(ceof)
  eof <- as.eof(ceof); attr(eof,"standard.error") <- NULL
  Z <- list(info='CORDEX runs',eof=eof,rcm.1=zoo(x1,order.by=index(ceof)))
  clim <- list(rcm.1=map.field(X,plot=FALSE))
  rcmnames <- attr(ceof,'model_id')
  gcmnames <- attr(ceof,'driving_model_id')
  gcmrip <- attr(ceof,'driving_model_ensemble_member')
  for (i in 1:attr(ceof,'n.apps')) {
    xi <- attr(ceof,paste('appendix.',i,sep=''))
    rcmnames <- c(rcmnames,attr(xi,'model_id'))
    gcmnames <- c(gcmnames,attr(xi,'driving_model_id'))
    gcmrip <- c(gcmrip,attr(xi,'driving_model_ensemble_member'))
    Z[[paste('rcm.',i+1,sep='')]] <- zoo(coredata(xi),order.by=index(xi))
    clim[[paste('rcm.',i+1,sep='')]] <- map.field(attr(X,paste('appendix.',i,sep="")))
  }
  attr(Z,'mean') <- clim
  attr(Z,'model_id') <- list(rcm=rcmnames,gcm=gcmnames,gcm_rip=gcmrip)
  class(Z) <- c('dsensemble','eof','list')
  ceof <- Z
  save(ceof,file=paste('ceof.rcm',varid,it,'rda',sep='.'))
  invisible(ceof)
}

as.field.commonEOFS <- function(x,verbose=FALSE) {
  if(verbose) print("as.field.ceof")
  Y <- as.field(x, anomaly=TRUE, verbose=verbose)
  Y0 <- Y
  for(i in 1:length(Y)) {
    clim.i <- coredata(attr(x,"mean")[[i]])
    Y.i <- coredata(Y[[i]]) 
    Y.i <- aperm(apply(Y.i,1,function(x) x + clim.i), c(2,1))
    coredata(Y[[i]]) <- Y.i
  }
  invisible(Y)
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
      #for(k in 1:(i-1)) Y.new[k,j] <- Y[k,]
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
    gcm.srex.name <- NA; gcm.srex.label <- NA
    
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
    
    if(!is.null(xx$model$srex)){
      for(mi in 1:length(xx$model$srex)){
        gcm.srex.label <- xx$model$srex[[mi]]$label
        gcm.srex.name <- xx$model$srex[[mi]]$name
      }
    }
    mx <- data.frame(project_id=project_id, url=url, filename=filename,
                     dim=paste(dim,collapse=","), dates=dates, var=paste(var,collapse=","),
                     longname=paste(longname,collapse=","), unit=paste(vunit,collapse=","),
                     #var_id=paste(vid,collapse=","), 
                     resolution=res, lon=lon.rng, lon_unit=lon.unit, lat=lat.rng, lat_unit=lat.unit,
                     experiment_id=experiment_id, frequency=frequency, 
                     creation_date=creation_date, #tracking_id=tracking_id,
                     gcm=gcm, gcm_rip=gcm.rip,
                     gcm.srex.label=gcm.srex.label, gcm.srex.name=gcm.srex.name)#, gcm_version=gcm.v, gcm_realm=gcm.realm)
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
