\name{commonEOF}
\alias{commonEOF}
\alias{commonEOF.gcm}
\alias{commonEOF.rcm}
\alias{subset.commonEOF}
\alias{map.commonEOF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate the common EOFs of data sets
}
\description{Computes common EOFs (a type of principal component
  analysis) for combinations of data sets. See \code{link[esd]{EOF.comb}}. 
  Here the method is applied to climate model results (GCMs and RCMs). 

  Select a subset of the common EOFs with \code{subset.commonEOF}.
  
  The common EOFs take less space than the original data. An estimate of 
  the original data sets (the ensemble mean or for an individual model) 
  can be calculated from the common EOFs with \code{link[map.commonEOF]}.
}
\usage{
  commonEOF.gcm(select=1:9,varid='tas',destfile=NULL,destfile.ceof=NULL,
                it='annual',is=NULL,verbose=FALSE)
  commonEOF.rcm(select=1:9,varid='tas',destfile=NULL,destfile.ceof=NULL,
                it='annual',is=NULL,verbose=FALSE)
  subset.commonEOF(x,it=NULL,is=NULL,ip=NULL,im=NULL,verbose=FALSE)
  map.commonEOF(x,it=NULL,is=NULL,ip=NULL,im=NULL,FUN="mean",plot=FALSE,
                colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,show=TRUE),
                verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{select}{index for selecting a subset of files/models.}
  \item{destfile}{The name of the local destination files for the original data.}
  \item{it}{A list or data.frame providing time index, e.g. year or month.}
  \item{is}{A list or data.frame providing space index in terms of lon and lat.}
  \item{ip}{A vector (integers) providing pattern index, i.e., which EOF patterns to use.}
  \item{im}{A vector (integers) providing model index, i.e., which of the models to use.}
  \item{varid}{variable name - see \code{link[esd]{retrieve}}}
  \item{verbose}{Logical value defaulting to FALSE. If FALSE, do not
    display comments (silent mode). If TRUE, displays extra information on
    progress.}
  \item{destfile.ceof}{The name of the local destination file of the commonEOF object.}
  \item{x}{A common EOF object, output of \code{commonEOF.rcm} or \code{commonEOF.gcm}.}
  \item{colbar}{The colour scales defined through \code{link[esd]{colscal}}. Users can
                specify the colour `pal'*ette (`pal'), the number of
                breaks (`n'), values of `breaks'. The `rev' argument,
                will produce a reversed color bar if set to TRUE.}
}
\value{A "dsensemble" "field" "list" object containing the common EOFs. The common EOF 
       is stored as a faux "dsensemble" object to take less space and so that we can use 
       some of the functions in esd designed for this type of object.}
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
  ## Calculate the common EOF of the first three CMIP5 GCMs.
  ## for the winter season (djf = december, january, february)
  ceof <- commonEOF.gcm(select=1:3,is=list(lon=c(-60,60),lat=c(45,80)),it="djf")
  
  ## Reproduce and plot a map of the original GCM data for the first model.
  X1 <- map.commonEOF(ceof,im=1)
  map(X1)
  
  ## Subset data from 1971-2000 and calculate the ensemble median for the period.
  ceof.7100 <- subset(ceof,it=c(1971,2000))
  X.median <- map.commonEOF(ceof.7100,FUN="median")
  map(X.median)
}
\keyword{ cds }
\keyword{ utilities }