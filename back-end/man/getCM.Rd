\name{getCM}
\name{cmip5.urls}
\name{cordex.urls}
\alias{getCM}
\alias{getGCMs}
\alias{getRCMs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCM(url=NULL,destfile='CM.nc',lon=NULL,lat=NULL,force=FALSE,verbose=FALSE)
cmip5.urls(experiment='rcp45',varid='tas',
           url="http://climexp.knmi.nl/CMIP5/monthly/", 
           path=NULL,off=FALSE,force=FALSE,verbose=FALSE)
cordex.urls(experiment='rcp45',varid='tas',
           url="http://climexp.knmi.nl/CMIP5/monthly/", 
           path=NULL,off=FALSE,force=FALSE,verbose=FALSE)
getatt(fname)
getGCMs(select=1:3,varid='tas',destfile=destfile=paste(rep('GCM',3),1:3,'.nc',sep=''))
getRCMs(select=1:3,varid='tas',destfile=destfile=paste(rep('RCM',3),1:3,'.nc',sep=''))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{url}{The URL locating the data files over Internet
%%     ~~Describe \code{url} here~~
}
  \item{destfile}{The name of the local destination file whewn downloading
%%     ~~Describe \code{destfile} here~~
}
  \item{lon}{longitude range to retrieve
%%     ~~Describe \code{lon} here~~
}
  \item{lat}{latitude range to retrieve
%%     ~~Describe \code{lat} here~~
}
  \item{fname}{file name}
  \item{select}{index for selecting a subset}
  \item{varid}{variable name - see \code{link[esd]{retrieve}}}
  \item{path}{local file path used for a speedier testing}
  \item{it}{Time index for subsetting - see \code{link[esd]{subset}}}
  \item{is}{Space index for subsetting - see \code{link[esd]{subset}}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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

## The function is currently defined as
function (url = "https://climexp.knmi.nl/CMIP5/monthly/tas/tas_Amon_ACCESS1-0_historical_000.nc", 
    destfile = "CM.nc", lon = NULL, lat = NULL) 
{
    lok <- try(download.file(url = url, destfile = destfile))
    if (inherits(loc, "try-error")) 
        return()
    X <- retrieve(destfile, lon = lon, lat = lat)
    cid <- getatt(destfile)
    cid$area.mean <- aggregate.area(X, FUN = "mean")
    cid$url <- url
    return(cid)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
