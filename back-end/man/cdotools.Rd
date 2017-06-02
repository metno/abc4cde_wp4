\name{cdo}
\alias(cdo.mean)
\alias(cdo.timeSd)
\alias(cdo.spatSd)
\alias(cdo.gridcor)
\alias(cdo.command)
\alias(gunzip)

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  CDO calculations in R
}
\description{Perform operations on netCDF files using CDO. More efficient than doing the same calculations directly in R.
}
\usage{
  cdo.mean(model.file,period=c(1981,2010),mask=NULL,seasonal=FALSE)
  cdo.timeSd(model.file,period=c(1981,2010),mask=NULL,seasonal=FALSE)
  cdo.spatSd(model.file,period=c(1981,2010),mask=NULL,seasonal=FALSE)
  cdo.gridcor(model.file,reference.file,period=c(1981,2010),mask=NULL,seasonal=FALSE)
  cdo.command(commands,input,infile,outfile,intern=FALSE)
  gunzip(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model.file}{Path to a netCDF file.}
  \item{period}{A list providing a time index in the form of a year range.}
  \item{mask}{A mask to apply to the data. If NULL use no mask. What can mask be except for NULL?}
  \item{seasonal}{TRUE: Calcuate seasonal values. FALSE: Caculate annual values.}
}
\value{Don't know what type of object the output is. The results are also saved in a netCDF file.}
\references{
  %% ~put references to the literature/web site here ~
}
\author{
  Kajsa Parding <kajsamp@met.no>, Olle Raety
}
\note{
  %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  \seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
  }
\examples{
}
\keyword{ cdo }
\keyword{ utilities }