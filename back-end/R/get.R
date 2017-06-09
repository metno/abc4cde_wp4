#Helper function to find files linux environment (add Windows counterpart)
find.file <- function(filename) {
  command <- paste("find -L $HOME -name",filename,sep=" ")
  fullpath <- system(command,intern=TRUE)
  if(length(fullpath)==0) return(FALSE)
  return(fullpath)
}

#Search and read a shapefile
get.shapefile <- function(filename,with.path=FALSE){
  fullname <- filename
  if(!with.path){
    fullname <- find.file(filename)
  }
  readOGR(fullname,verbose=FALSE)
}

# Provides path to local files with reference data
getReference <- function(reference,variable,verbose=FALSE) {
  if(verbose) print("getReference")
  #path <- system("echo $EXTERNAL_DATA",intern=TRUE)
  file.name <- switch(paste(reference,variable,sep="."),
                      era.tas="era-interim_monthly_1979-2016_tas.2.5deg.nc",
                      era.pr="era-interim_monthly_1979-2016_pr.2.5deg.nc",
                      cfsr.tas="cfsr_tmp2m_mon.nc",
                      cfsr.pr="cfsr_prate_mon.nc",
                      eobs.tas="tg_0.50deg_reg_v14.0_mon.nc",
                      eobs.pr="rr_0.50deg_reg_v14.0_mon.nc")
  file.name <- find.file(file.name)[1]
  invisible(file.name)
  #invisible(paste(path,file.name,sep="/"))
}

# Generate name of GCM or RCM file path/prefix.number.variable.nc
get.name <- function(number,variable,is.rcm=FALSE,verbose=FALSE) {
  if(verbose) print("get.name")
  path <- system("echo $EXTERNAL_DATA",intern=TRUE)
  prefix <- "GCM"
  if(is.rcm) prefix <- "RCM"
  file.name <- paste(paste(prefix,number,sep=""),variable,"nc",sep=".")
  invisible(paste(path,file.name,sep="/"))
}

# Call to a python script which downloads data from the public ECMWF data server
python.getEra <- function(start,end,variable,steps,type,stream,outfile,verbose=FALSE) {
  if(verbose) print("python.getEra")
  #script <- "python python/getMonthlyERA.py"
  path.era <- find.file("getMonthlyERA.py")
  script <- paste("python",path.era[1])
  if(verbose) print("Run ECMWF Python script for downloading ERA data")
  system.command <- paste(script," -f ",start," -l ",end," -v ",variable,
                          " -s ",steps," -t ",type," -r ",stream," -o ",outfile, sep="")
  system(system.command,wait=TRUE)
  if(verbose) print("--- end python.getEra")
}

#Retrieve monthly data for 2m temperature and precipitation from the ECMWF public repository.
#The use of this function requires that the ECMWF key and python libraries are installed on the machine.
#See instructions in https://software.ecmwf.int/wiki/display/WEBAPI/Access+ECMWF+Public+Datasets
#The function also requires that cdo is installed on the operating computer.
getERA <- function(variable.name,start=1979,end=2016,griddes="cmip_1.25deg_to_2.5deg.txt",
                   destfile=NULL,verbose=FALSE){
  if(verbose) print("getERA")
  griddes <- find.file(griddes)[1]
  if(any(match(c("tas","tmp","temp","temperature","t2m"),variable.name,nomatch=0))) {
    if(verbose) print("variable: temperature")
    varID<- "167.128"
    stream <- "moda"
    steps <- "0"
    type <- "an"
    commands <- c("-f","nc","-copy","-remapcon","-chname")
    input <- c("","","",griddes,"2t,tas")
  } else if(any(match(c("pre","prc","prec","precipitation","pr"),variable.name,nomatch=0))) {
    if(verbose) print("variable: precipitation")
    varID<- "228.128"
    stream <- "mdfa"
    steps <- "12-24/24-36" # Select the 24 and 36h 12-hour long forecast in order to reduce the spin-up effect.
    type <- "fc"
    commands <- c("-f","nc","-copy","-monsum","-remapcon","-chname")
    input <- c("","","","",griddes,"2t,tas")
  }
  #griddes <- find.file(griddes)
  if(is.null(destfile)) destfile <- paste("era-interim_monthly_",paste(start,end,sep="-"),"_",variable.name,".grib",sep="")
  if(!file.exists(destfile)) {
    if(verbose) print("File does not exist. Download with ECMWF Python tool.")
    python.getEra(start, end, varID, steps, type, stream, destfile, verbose=verbose)
  }
  outfile <- paste(gsub('.{5}$', '',destfile),"2.5deg",'nc',sep=".")
  if(!file.exists(outfile)) {
    if(verbose) print("File with 2.5deg data does not exist. Regrid with CDO.")
    cdo.command(commands,input,infile=destfile,outfile=outfile)
  }
  if(verbose) print("Retrieve data from netCDF file.")
  X <- retrieve(outfile)
  cid <- getatt(outfile)
  if(verbose) print("Calculate area mean and sd.")
  cid$area.mean <- aggregate.area(X,FUN='mean')
  cid$area.sd <- aggregate.area(X,FUN='sd')
  cid$url <- NA
  cid$dates <- paste(range(index(X)),collapse=",")
  if(verbose) print("Get information about the model and netCDF file.")
  ncid <- nc_open(outfile)
  model <- ncatt_get(ncid,0)
  nc_close(ncid)
  cid$model <- model
  cid$project_id <- cid$model$project_id
  cid$srex <- get.srex.region(outfile,region=NULL,print.srex=FALSE,verbose=FALSE)
  if(verbose) print("--- end getERA ---")
  return(cid)
}

#A test function to retrieve CRU data from CEDA databases.
getCRU <- function(username,passwd,variable="tmp",version="4.00",
                   griddes="cmip_1.25deg_to_2.5deg.txt",destfile=NULL,
                   time.resol=NULL, verbose=FALSE) {
  if(verbose) print("getCRU")
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
    commands <- c("-f","nc","-copy","-remapcon")
    input <- c("","","",griddes)
    cdo.command(commands,input,infile=destfile,outfile=outfile) 
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
getCFSR <- function(variable="tas",destfile=NULL,lon=NULL,lat=NULL,
                    griddes="cmip_1.25deg_to_2.5deg.txt",verbose=FALSE) {
  if(verbose) print("getCFSR")
  url.path <- "https://climexp.knmi.nl/CFSR"
  griddes <- find.file(griddes)
  if(variable=="tas"){
    filename <- "cfsr_tmp2m.nc"
    commands <- c("-f","nc","-copy","-remapcon","-monavg","-chname")
    input <- c("","","",griddes,"","TMP_2maboveground,tas")
  } else if(variable=="pr"){
    filename <- "cfsr_prate.nc"
    commands <- c("-f","nc","-copy","-remapcon","-monavg","-chname")
    input <- c("","","",griddes,"","PRATE_surface,pr")
  }
  if(!file.exists(filename)) download.file(paste(url.path,filename,sep="/"),destfile=filename)
  if(is.null(destfile)) destfile <- paste(sub("\\.[[:alnum:]]+$", "", filename, perl=TRUE),"mon.nc",sep="_")
  if(!file.exists(destfile)) cdo.command(commands,input,infile=filename,outfile=destfile)
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
  return(cid)
}

#Get daily EOBS data and convert it to monthly averages. Version and resolution
#selection not implemented yet.
getEOBS <- function(variable="tas", destfile=NULL, resolution="0.50", version="14",
                    verbose=FALSE) {
  if(verbose) print("getEOBS")
  url.path <- "http://www.ecad.eu/download/ensembles/data/Grid_0.50deg_reg"
  if(variable=="tas"){
    filename <- "tg_0.50deg_reg_v14.0.nc.gz"
  }else if(variable=="pr"){
    filename <- "rr_0.50deg_reg_v14.0.nc.gz"
  }else{
    return("Not implemented yet!")
  }
  if(!file.exists(filename)) download.file(paste(url.path,filename,sep="/"),destfile=filename)
  gunzip(filename)
  filename <- sub("\\.[[:alnum:]]+$", "", filename, perl=TRUE)
  if(is.null(destfile)) destfile <- paste(paste(system("echo $EXTERNAL_DATA",intern=T),sub("\\.[[:alnum:]]+$", "", filename, perl=TRUE),sep="/"),"mon.nc",sep="_")
  print(destfile)
  commands <- c("-f","nc","-copy","-monavg")
  input <- c("","","","")
  if(!file.exists(destfile)) cdo.command(commands,input,infile=filename,outfile=destfile)
  X <- retrieve(destfile,lon=lon,lat=lat,verbose=verbose)
  cid <- getatt(destfile) 
  cid$url <- paste(url.path,filename,sep="/")
  cid$area.mean <- aggregate.area(X,FUN='mean',na.rm=T)
  cid$area.sd <- aggregate.area(X,FUN='sd',na.rm=T)
  ncid <- nc_open(destfile)
  model <- ncatt_get(ncid,0)
  nc_close(ncid)
  cid$model <- model
  return(cid)
}
