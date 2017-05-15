## ABC4CDE/DEMC - R-script for prototype tool WP4 Evaluation and Quality Control (EQC)
## Rasmus Benestad. Blindern, Oslo 2016-10-17

## A cloolection of R-scripts which will be used in a prototype tool that runs regularly in the background
## and collects metadata and performs evaluation of the data. 

## Evaluation of CMIP ensemble based on NCEP/NCAR 1, ERAINT, and MERRA reanalyses
## Surface temperature. Rank-test for grid-points values

testrank <- function(x,m,N,d1) {
  ## Find the rank number of m observations
  ## Need to unwrap the dimensions
  dim(x) <- c(N,d1)
  number <- apply(x,1,function(x) order(x)[1:m])
  return(c(number))
}


EQC.gcmbias <- function(x,...) {

}

## This function evaluates whole ensembles and compares the sample of simulated results to corresponding
## 'observaitons' (one or several reanalyses). The testing is done through rank-statistics, and we expect
## thatt the rank number follows a uniform distribution if the simulated results and observations blong to
## the same statistical population. 
ECQ.ensemble <- function(obs=c('air.mon.mean.nc','ETAINT_t2m.mon.nc','MERRA'),
                         path='CMIP5.monthly',pattern='tas_',it=c(1980,2015),is=NULL,
                         anomaly=FALSE) {
  require(esd)
  
  ## Get all the GCM data
  gcms <- list.files(path=path,pattern=pattern,full.names=TRUE)
  n <- length(gcms)
  N <- n + length(obs)
  m <- length(obs)
  ## Use the grid of the first reanalysis and interploate the other fields onto this grid
  y1 <- subset(retrieve(obs[1]),it=it,is=is)
  if (anomaly) y1 <- anomaly(y1)
  d <- dim(y1)
  
  ## Set up a data matrix containing the data - start with the first reanalysis
  X <- rep(NA,N*d[1]*d[2]); dim(X) <- c(N,d)
  cnames <- rep("",N)
  X[1,,] <- coredata(y1)
  cnames[1] <- paste(attr(y1,'model_id'),attr(y1,'run'))
  
  if (m > 1) {
    for (i in 2:m) {
      y <- subset(retrieve(obs[i]),it=it,is=is)
      y <- regrid(y,is=y1)
      if (anomaly) y <- anomaly(y)
      X[i,,] <- coredata(y)
      cnames[i] <- paste(attr(y,'model_id'),attr(y,'run'))
    }
  }
  
  ## Do the GCMs:
  for (i in 1:n) {
    y <- subset(retrieve(gcms[i]),it=it,is=is)
    y <- regrid(y,is=y1)
    if (anomaly) y <- anomaly(y)
    X[i+m,,] <- coredata(y)
    cnames[i+m] <- paste(attr(y,'model_id'),attr(y,'run'))
  }
  
  ## convert the matrix to 2D
  dim(X) <- c(N*d[1],d[2])
  
  ## For each gridpoint, accumulate the ranks for all the time steps.
  ## The dimension of the results should be [m*d[1],d[2]], where te elements of the 
  ## first dimension is expected to follow a uniform distribution if the GCMs and reanalysis
  ## belong to the same population.
  testmap <- apply(X,1,testrank,m,N,d[1])
}


EQC.stations <- function() {
  ## wilcox.test for assessment against station values
} 

EQC.commonEOF <- function() {
  ## For visualisation and assessment/comparison of spatio-temporal covariance structure
}

EQC.ensemblenorm <- function() {
  ## Compare the grid points to a normal distribution:  the Kolmogorov–Smirnov
  ## https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
}

compare.fields <- function(x,y=NULL,lplot=FALSE,type=c("correlation","rmsd"),
                           filename=NULL,verbose=FALSE,...) {
  if(verbose) print("compare.fields")
  if(!is.null(y) & inherits(x,"field") & inherits(y,"field")) {
    if(verbose) print("Calculate comparative statistics")
    x <- subset(x,is=list(lon=range(lon(y)),lat=range(lat(y))))
    y <- subset(y,is=list(lon=range(lon(x)),lat=range(lat(x))))
    if(inherits(x,"annual") & !inherits(y,"annual")) y <- annual(y)
    if(inherits(y,"annual") & !inherits(x,"annual")) x <- annual(x)
    if(inherits(x,"seasonal") & !inherits(y,"seasonal")) y <- subset(as.4seasons(y),it=season(x)[1])
    if(inherits(y,"seasonal") & !inherits(x,"seasonal")) x <- subset(as.4seasons(x),it=season(y)[1])
    x <- subset(x,it=index(x) %in% index(y))
    y <- subset(y,it=index(y) %in% index(x))
    if(verbose) print("Calculate correlation")
    r <- corfield(x,y,plot=FALSE)
    attr(r,"variable") <- "correlation"
    attr(r,"longname") <- "correlation between fields"
    attr(r,"unit") <- "-1 to -1"
    if(verbose) print("Calculate difference in means")
    mdiff <- apply(x,2,mean)-apply(y,2,mean)
    mdiff <- attrcp(r,mdiff)
    attr(mdiff,"variable") <- "meandiff"
    attr(mdiff,"longname") <- "difference in means"
    attr(mdiff,"unit") <- attr(x,"unit")
    class(mdiff) <- "corfield"
    if(verbose) print("Calculate difference in trend")
    fn <- function(x) if(any(!is.na(x))) return(trend.coef(x)) else return(NA)
    tdiff <- apply(x,2,fn) - apply(y,2,fn)
    tdiff <- attrcp(r,tdiff)
    attr(tdiff,"variable") <- "trenddiff"
    attr(tdiff,"longname") <- "difference in trend"
    attr(tdiff,"unit") <- paste(attr(x,"unit"),"/decade",sep="")
    class(tdiff) <- "corfield"
    if(verbose) print("Calculate RMSD and normalised RMSD")
    rmsd <- apply((x-y)^2,2,function(x) sqrt(sum(x,na.rm=TRUE)/sum(!is.na(x))))
    nrmsd <- 100*rmsd/apply(x,2,function(x) diff(range(x,na.rm=TRUE)))
    rmsd <- attrcp(r,rmsd)
    attr(rmsd,"variable") <- "RMSD"
    attr(rmsd,"longname") <- "root mean square deviation"
    attr(rmsd,"unit") <- attr(x,"unit")
    class(rmsd) <- "corfield"
    nrmsd <- attrcp(r,nrmsd)
    attr(nrmsd,"variable") <- "NRMSD"
    attr(rmsd,"longname") <- "root mean square deviation normalised by the range of the data"
    attr(nrmsd,"unit") <- "%"
    class(nrmsd) <- "corfield"
    z <- list(correlation=r,meandiff=mdiff,trenddiff=tdiff,rmsd=rmsd,nrmsd=nrmsd)
  } else {
    z <- x
  }
  if(lplot) {
    if(verbose) print("Plot comparison between fields")
    stopifnot(inherits(z,"corfield") |
             (inherits(z,"list") & inherits(z[[1]],"corfield")))
    if(inherits(z,"corfield")) {
      eval(parse(text=paste("z <- list(",attr(z,"variable"),"=z)",sep="")))
      type <- c(attr(z,"variable"))
    }
    type <- type[type %in% names(z)]
    if(is.null(type)) type <- names(z)
    if(!is.null(filename)) {
      pdf(filename, 3.2*length(type), 4.0)
    } else {
      dev.new(width=3.2*length(type),height=4.0)
    }
    par(mar=c(4.5,2.5,3.5,0.5),mgp=c(1.5,0.5,0))
    for (i in 1:length(type)) {
      z.i <- z[[which(names(z)==type[[i]])]]
      cb.i <- select.colbar(z.i)
      fig.i <- c((i-1)/length(type),i/length(type),0,1)
      if(i==1) par(fig=fig.i) else par(fig=fig.i,new=TRUE)
      map(z.i,colbar=cb.i)
    }
    if(!is.null(filename)) dev.off()
  }
  invisible(z)
}
