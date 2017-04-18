## ABC4CDE/DECM WP4 R-scripts for propotype climate data storage (CDS) tool for 
## processing data, metadata and evaluationand quality control (EQC)
## Rasmus.Benestad@met.no 2017-02-15

## Function that extracts metadata stored in lists and sublists so that it's easier to search and
## organise into tables/spreadsheets. 
##
##metaxtract <- function(x=NULL) {
##  ## argument 'x' is input from getGCMs, getRCMs, testGCM, etc
## if (is.null(x)) x <- c(getGCMs(),getRCMs())
##  
##  gcms <- names(x); n <- length(gcms)
##  print(gcms)
##  
##  for (i in 1:n) {
##    xx <- x[[gcms[i]]]                       ## Tidier - use xx for the next level in the list structure
##    meta <- names(xx); m <- length(meta)
##    if (i==1) {
##      X <- matrix(rep("NA",n*m),n,m) ## set up a matrix
##      colnames(X) <- meta; rownames(X) <- gcms
##    }
##    for (ii in 1:m) {
##      if (is.list(xx[[ii]])) y <- unlist(xx[[ii]]) else 
##      if (inherits(xx[[ii]],'station')) y <- as.character(round(mean(subset(xx[[ii]],it=c(2070,2100))),2)) else
##      y <- xx[[ii]] 
##      print(class(xx[[ii]])); print(y)
##      if (length(y)>1) y <- paste(y,collapse='-')
##      X[i,ii] <- y
##    }
##  }
##  X -> metaxtract
##  save(metaxtract,file='metaxtract.rda')
##  return(X)
##}
