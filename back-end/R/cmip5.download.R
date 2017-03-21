# Description / download.cmip5 !

# This script allows downloading CMIP5 climate variables automatically from 
#  the climate explorer KNMI Web link "http://climexp.knmi.nl"

## Created by           abdelkader.mezghani@met.no 
## Created              2012-09-18
## Updated by           kajsa.parding@met.no
## last updates         2016-09-27
## revised              rasmus.benestad@met.no 2017-02-14

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
