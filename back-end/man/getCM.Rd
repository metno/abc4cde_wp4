\name{getCM}
\alias{getCM}
\alias{getGCMs}
\alias{getRCMs}
\name{cmip5.urls}
\name{cordex.urls}
\name{thredds.urls}
\alias{getatt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Download GCM and RCM data
}
\description{\code{getCM} downloads GCM or RCM data from given a URL. \code{getGCMs} and \code{getRCMs} retrieves CMIP5 and CORDEX data from the KNMI climate data explorer, generating the URLs with \code{cmip5.urls} and \code{cordex.urls}. \code{getatt} reads and extracts the attribute information in a netCDF files and stores this in a list object.
}
\usage{
getCM(url=NULL,destfile='CM.nc',lon=NULL,lat=NULL,force=FALSE,verbose=FALSE)
getatt(url=NULL,destfile='CM.nc',lon=NULL,lat=NULL,force=FALSE,verbose=FALSE)
cmip5.urls(experiment='rcp45',varid='tas',
           url="http://climexp.knmi.nl/CMIP5/monthly/", 
           path=NULL,off=FALSE,force=FALSE,verbose=FALSE)
cordex.urls(experiment='rcp45',varid='tas',
           url="http://climexp.knmi.nl/CMIP5/monthly/", 
           path=NULL,off=FALSE,force=FALSE,verbose=FALSE)
thredds.urls(url.rel="raw/tas",pattern=".*EUR-11.*.nc",select=NULL,
           url.catalog="http://thredds.met.no/thredds/catalog/postclim/data/CORDEX-EUR11",
           url.data="http://thredds.met.no/thredds/dodsC/postclim/data/CORDEX-EUR11",
           verbose=FALSE)
getGCMs(select=1:3,varid='tas',destfile=destfile=paste(rep('GCM',3),1:3,'.nc',sep=''))
getRCMs(select=1:3,varid='tas',destfile=destfile=paste(rep('RCM',3),1:3,'.nc',sep=''))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{url}{The URL locating the data files over Internet}
  \item{destfile}{The name of the local destination file when downloading}
  \item{fname}{The filename of a local netCDF file}
  \item{lon}{longitude range to retrieve}
  \item{lat}{latitude range to retrieve}
  \item{select}{index for selecting a subset of files/models to download}
  \item{varid}{variable name - see \code{link[esd]{retrieve}}}
  \item{path}{local file path to the destination file}
  \item{it}{Time index for subsetting - see \code{link[esd]{subset}}}
  \item{is}{Space index for subsetting - see \code{link[esd]{subset}}}
  \item{url.rel}{URL relative to \code{url.base} (and \code{url.download)}}
  \item{url.base}{Base URL to a catalog of available files.}
  \item{url.download}{Base URL to location where data can be downloaded. If NULL, url.base is used. (I had problems listing and downloading files from the same location on the metno THREDDS server and solved it by including two different base URLs.)}
  \item{force}{Logical value defaulting to FALSE. If FALSE, do not
    download data from source if a local file \code{path/destfile} is available. If TRUE, download data and overwrite local file.}
  \item{verbose}{Logical value defaulting to FALSE. If FALSE, do not
    display comments (silent mode). If TRUE, displays extra information on
    progress.}
}
\value{
\code{getCM} downloads the data and returns a "zoo" "field" object. 
\code{getGCMs} and \code{getRCMs} returns a list of "zoo" "field" objects corresponding to different models.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Kajsa Parding <kajsamp@met.no>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

x <- getCM(url='https://climexp.knmi.nl/CMIP5/monthly/tas/tas_Amon_ACCESS1-0_historical_000.nc',
      destfile='GCM1.nc',lon=c(-60,60),lat=c(45,80),force=TRUE,verbose=TRUE)
    
x.info <- getatt(x)
y <- getGCMs(select=1:3,lon=c(-60,60),lat=c(45,80))
map(y[[1]])
}

\keyword{ cds }
\keyword{ utilities }