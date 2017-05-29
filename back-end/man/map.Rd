\name{map.ensemble}
\alias{map.ensemble}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Map climate model projections.
}
\description{
Produce maps showing the average or change in climate based on climate model projections.
}
\usage{
map.ensemble(ceof,im=NULL,ip=NULL,is=NULL,type=NULL,new=TRUE,FUN="mean",
             colbar=list(pal=NULL,breaks=NULL,show=TRUE,rev=FALSE),
             verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ceof}{A common EOF object, output of \code{commonEOF.rcm} or \code{commonEOF.gcm}.}
  \item{is}{A list or data.frame providing space index in terms of lon and lat.}
  \item{ip}{A vector (integers) providing pattern index, i.e., which EOF patterns to use.}
  \item{im}{A vector (integers) providing model index, i.e., which of the models to use.}
  \item{type}{Type of map to produce. The valid options are: 
              "ccff" and "ccnf" (climate change far future and near future), 
              "mvff", "mvnf", "mvpd" (mean value far future, near future and present day).}
  \item{verbose}{Logical value defaulting to FALSE. If FALSE, do not
    display comments (silent mode). If TRUE, displays extra information on
    progress.}
  \item{colbar}{The colour scales defined through \code{link[esd]{colscal}}. Users can
                specify the colour `pal'*ette (`pal'), the number of
                breaks (`n'), values of `breaks'. The `rev' argument,
                will produce a reversed color bar if set to TRUE.}
%%     ~~Describe \code{x} here~~
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{x}{A "field" "zoo" object containing the average climate or climate change depending on the selected type.}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ visualisation }% use one of  RShowDoc("KEYWORDS")