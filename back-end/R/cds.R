## ABC4CDE/DEMC - R-script for prototype tool WP4
## Rasmus.Benestad@met.no  Oslo, Norway, 2017-02-14
##

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
