\name{metaextract}
\alias{metaextract}
\alias{metaextract.cmip5}
\alias{metaextract.cordex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract metadata
}
\description{
Extract metadata for GCMs and RCMs.
}
\usage{
metaextract(x)
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
    gcms <- names(x)
    n <- length(gcms)
    for (i in 1:n) {
        xx <- x[[gcms[i]]]
        meta <- names(xx)
        m <- length(meta)
        if (i == 1) {
            X <- matrix(rep("NA", n * m), n, m)
            colnames(X) <- meta
            rownames(X) <- gcms
        }
        for (ii in 1:m) {
            if (is.list(xx[[ii]])) 
                y <- unlist(xx[[ii]])
            else if (inherits(xx[[ii]], "station")) 
                y <- as.character(round(mean(subset(xx[[ii]], 
                  it = c(2070, 2100))), 2))
            else y <- xx[[ii]]
            print(class(xx[[ii]]))
            print(y)
            if (length(y) > 1) 
                y <- paste(y, collapse = "-")
            X[i, ii] <- y
        }
    }
    return(X)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
