## ABC4CDE/DEMC - R-script for prototype tool WP4
## andreas.dobler@met.no  Oslo, Norway, 2017-05-16
##

##  Function of an R-package that retrieves metadata from files located on OPeNDAP servers.
require(PCICt) ##Used to get the time axes right for 360 day and non-leapyear calendars.

#Function to extract the metadata form NetCDF files on OPeNDAP servers
metaextract_opendap <- function(dataurl=NULL, verbose=FALSE) {
  if(verbose) print("metaextract_opendap")
  
  #Make connection to file and get global attributes, variable and dimension lists
  nc <- nc_open(dataurl)
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
  return(list(globat=globat,varlist=varlist,dimlist=dimlist,dataurl=dataurl,range=range(timeline)))
}
