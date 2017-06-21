## Compute the common EOFs for GCMs and save the results for the front-end
commonEOF.gcm <- function(select=1:9,varid='tas',destfile=NULL,destfile.ceof=NULL,
                          it='annual',is=NULL,verbose=FALSE) {
  if(verbose) print("commonEOF.gcm")
  if(is.null(destfile)) destfile <- paste('GCM',select,'.',varid,'.nc',sep='')
  getGCMs(select=select,varid=varid,destfile=destfile)
  X <- NULL
  for (fname in destfile) {
    if(verbose) print(paste("retrieve",fname))
    x <- retrieve(fname,verbose=verbose)
    if (!is.null(it)) {
      if (tolower(it)=='annual') x <- annual(subset(x,is=is),verbose=verbose) else
        x <- subset(x,it=it,is=is,verbose=verbose)
    }
    if (is.null(X)) X <- x else X <- combine(X,x,verbose=verbose)
  }
  if(verbose) print("Calculate common EOF")
  ceof <- EOF(X,verbose=verbose)
  plot(ceof)
  
  ## Need to reformat the ceof-object to fit the set-up for the R-shiny app in the front-end.
  if(verbose) print("Reformat the common EOF object")
  x1 <- coredata(ceof); attributes(x1) <- NULL; dim(x1) <- dim(ceof)
  eof <- as.eof(ceof); attr(eof,"standard.error") <- NULL
  Z <- list(info='CMIP5 runs',eof=eof,rcm.1=zoo(x1,order.by=index(ceof)))
  clim <- list(rcm.1=map.field(X,plot=FALSE))
  gcmnames <- attr(ceof,'model_id')
  gcmrip <- attr(ceof,'parent_experiment_rip')
  for (i in 1:attr(ceof,'n.apps')) {
    x1 <- attr(ceof,paste('appendix.',i,sep=''))
    Z[[paste('rcm.',i+1,sep='')]] <- zoo(coredata(x1),order.by=index(x1))
    gcmnames <- c(gcmnames,attr(x1,'model_id'))
    gcmrip <- c(gcmrip,attr(x1,'parent_experiment_rip'))
    clim[[paste('rcm.',i+1,sep='')]] <- map.field(attr(X,paste('appendix.',i,sep='')),plot=FALSE)
  }
  attr(Z,'mean') <- clim
  attr(Z,'model_id') <- list(gcm=gcmnames,gcm_rip=gcmrip)
  class(Z) <- c('dsensemble','eof','list')
  ceof <- Z
  if(is.null(destfile.ceof)) {
    destfile.ceof <- paste('ceof.gcm',varid,it,sep=".")
    if(!is.null(is)) fname <- paste(destfile.ceof,".",paste(round(is$lon),collapse="-"),"E.",
                                    paste(round(is$lat),collapse="-"),"N",sep="")
    destfile.ceof <- paste(destfile.ceof,"rda",sep='.')
  }
  save(ceof,file=destfile.ceof)
  invisible(ceof)
}

## Compute the common EOFs for RCMs save the results for the front-end
commonEOF.rcm <- function(select=1:9,varid='tas',destfile=NULL,destfile.ceof=NULL,
                          it='annual',is=NULL,verbose=FALSE) {
  if(verbose) print("commonEOF.rcm")
  if(is.null(destfile)) destfile <- paste(rep('CM',length(select)),select,'.',varid,'.nc',sep='')
  getRCMs(select=select,varid=varid,destfile=destfile,verbose=verbose)
  X <- NULL
  for (fname in destfile) {
    if(verbose) print(paste("retrieve",fname))
    x <- retrieve(fname)
    if (!is.null(it)) {
      if (tolower(it)=='annual') x <- annual(subset(x,is=is),verbose=verbose) else
        x <- subset(x,it=it,is=is,verbose=verbose)
    }
    if (is.null(X)) X <- x else X <- combine(X,x,verbose=verbose)
  }
  if(verbose) print("Calculate common EOF")
  ceof <- EOF(X,verbose=verbose)
  plot(ceof)
  
  ## Need to reformat the ceof-object to fit the set-up for the R-shiny app in the front-end.
  if(verbose) print("Reformat the common EOF object")
  x1 <- coredata(ceof); attributes(x1) <- NULL; dim(x1) <- dim(ceof)
  eof <- as.eof(ceof); attr(eof,"standard.error") <- NULL
  Z <- list(info='CORDEX runs',eof=eof,rcm.1=zoo(x1,order.by=index(ceof)))
  clim <- list(rcm.1=map.field(X,plot=FALSE))
  rcmnames <- attr(ceof,'model_id')
  gcmnames <- attr(ceof,'driving_model_id')
  gcmrip <- attr(ceof,'driving_model_ensemble_member')
  for (i in 1:attr(ceof,'n.apps')) {
    xi <- attr(ceof,paste('appendix.',i,sep=''))
    rcmnames <- c(rcmnames,attr(xi,'model_id'))
    gcmnames <- c(gcmnames,attr(xi,'driving_model_id'))
    gcmrip <- c(gcmrip,attr(xi,'driving_model_ensemble_member'))
    Z[[paste('rcm.',i+1,sep='')]] <- zoo(coredata(xi),order.by=index(xi))
    clim[[paste('rcm.',i+1,sep='')]] <- map.field(attr(X,paste('appendix.',i,sep="")),plot=FALSE)
  }
  attr(Z,'mean') <- clim
  attr(Z,'model_id') <- list(rcm=rcmnames,gcm=gcmnames,gcm_rip=gcmrip)
  class(Z) <- c('dsensemble','eof','list')
  ceof <- Z
  if(is.null(destfile.ceof)) {
    destfile.ceof <- paste('ceof.rcm',varid,it,sep=".")
    if(!is.null(is)) destfile.ceof <- paste(destfile.ceof,".",paste(round(is$lon),collapse="-"),"E.",
                                            paste(round(is$lat),collapse="-"),"N",sep="")
    destfile.ceof <- paste(destfile.ceof,"rda",sep='.')
  }
  save(ceof,file=destfile.ceof)
  invisible(ceof)
}

subset.commonEOF <- function(x,it=NULL,is=NULL,ip=NULL,im=NULL,verbose=FALSE) {
  if(verbose) print("subset.commonEOF")
  Y <- subset.dsensemble.multi(x,it=it,is=is,ip=ip,im=im,verbose=verbose)
  if(is.null(im)) im <- seq(length(x)-2)
  Y.mean <- attr(x,"mean")[im]
  if(length(im)==1) Y.mean <- list(Y.mean)
  if(!is.null(is)) {
    ok.lon <- attr(x[[2]],"longitude")>=min(is$lon) & 
      attr(x[[2]],"longitude")<=max(is$lon)
    ok.lat <- attr(x[[2]],"latitude")>=min(is$lat) & 
      attr(x[[2]],"latitude")<=max(is$lat)
    for(i in seq(1,length(Y.mean))) {
      clim.i <- Y.mean[[i]][ok.lon,ok.lat]
      Y.mean[[i]] <- clim.i
    }
  }
  attr(Y,"mean") <- Y.mean
  return(Y)
}

map.commonEOF <- function(x,it=NULL,is=NULL,ip=NULL,im=NULL,FUN=NULL,plot=FALSE,
                          colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,show=TRUE),
                          verbose=FALSE) {
  if(verbose) print("map.commonEOF")
  x <- subset.commonEOF(x,it=it,is=is,ip=ip,im=im,verbose=verbose)
  Y <- map(x,it=it,anomaly=TRUE,plot=FALSE,FUN=FUN,verbose=verbose)
  if(is.null(FUN)) FUN <- "" ## If FUN = NULL the following line doesn't work:
  if( FUN %in% c("mean","median","q5","q95") ) {
    if(is.null(im)) im <- seq(length(x)-2)
    clim <- unlist(coredata(attr(x,"mean")))
    dim(clim) <- c(ncol(Y),length(im))
    clim.avg <- apply(clim,1,mean)
    Y <- Y + aperm(array(rep(clim.avg,nrow(Y)),rev(dim(Y))),c(2,1))
  }
  if(plot) map(Y)
  invisible(Y)
}
