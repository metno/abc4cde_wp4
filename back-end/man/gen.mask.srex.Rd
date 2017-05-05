\name{gen.mask.srex}
\alias{getCM}
\alias{getatt}
\alias{getGCMs}
\alias{testGCM}
\alias{getRCMs}
\alias{cmip5.urls}
\alias{commonEOFS.gcm}
\alias{commonEOFS.rcm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gen.mask.srex(infile, region, ind = FALSE, inverse = FALSE, mask.values = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{Netcdf file from which the raster is created.
%%     ~~Describe \code{infile} here~~
}
  \item{region}{SREX region, for which the mask will be created
%%     ~~Describe \code{region} here~~
}
  \item{ind}{logical indicating whether only the indices of the pixels within the selected region are returned
%%     ~~Describe \code{ind} here~~
}
  \item{inverse}{if TRUE, the function returns the inverse mask for the selected SREX region as a raster object or the corresponding pixel indices, if ind = TRUE.
%%     ~~Describe \code{inverse} here~~
}
  \item{mask.values}{Which values are given for the mask. The mask values are set to 1 by default.
%%     ~~Describe \code{mask.values} here~~
}
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
#Not run:
tas.cmip <- getGCMs(1,"tas")
mask <- gen.mask.srex(infile="GCM1.tas.nc",region="ALA")
plot(mask)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
