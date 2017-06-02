\name{EQC}
\alias{EQC.gcmbias}
\alias{EQC.ensemble}
\alias{EQC.stations}
\alias{EQC.commonEOF}
\alias{EQC.ensemblenorm}
\alias{EQC.scatterplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluation and quality control
}
\description{
Various functions for evaluation and quality control of climate model data.

\code{EQC.gcmbias} calculates the bias of GCMs with regards to what? Not done.

\code{EQC.ensemble} compares the range of a GCM ensemble to the natural variability seen in observations.

\code{EQC.stations} evaluates climate model data against corresponding surface based observations.

\code{EQC.commonEOF} is a function for visualisation and assessment/comparison of spatio-temporal covariance structure.

\code{EQC.ensemblenorm} compares the grid points to a normal distribution: the Kolmogorov–Smirnov

\code{EQC.scatterplot} shows the average change of temperature and precipiation from the present day to far future for all members of a GCM or RCM ensemble. 

}
\usage{
EQC.gcmbias(x)
EQC.ensemble(x)
EQC.stations(x)
EQC.commonEOF(x)
EQC.ensemblenorm(x)
EQC.scatterplot(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
Kajsa Parding, kajsamp@met.no
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
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
