
maptype <- function(type="Climate change far future") {
  if(!is.null(type)) {
    switch(type,"Climate change near future"="ccnf",
           "Climate change far future"="ccff",
           "Mean value present day"="mvpd",
           "Mean value near future"="mvnf",
           "Mean value far future"="mvff")
  } else {
    "ccff"
  }
}

map.ensemble <- function(ceof,im=NULL,ip=NULL,is=NULL,type=NULL,new=TRUE,FUN="mean",
                         colbar=list(pal=NULL,breaks=NULL,show=TRUE,rev=FALSE),
                         verbose=FALSE) {
  if(verbose) print("map.ensemble")
  mt <- maptype(type)
  if(grepl("cc",mt)) {
    if(grepl("t2m|tas|temp",attr(ceof[[2]],"variable"))) {
      colbar$breaks <- seq(-10,10,1)
      if(is.null(colbar$pal)) colbar$pal <- "burd"
    }
    if(grepl("pr",attr(ceof[[2]],"variable"))) {
      colbar$breaks <- seq(-2,2,0.2)
      if(is.null(colbar$pal)) {colbar$pal <- "burd"; colbar$rev <- TRUE} 
    }
    it1 <- c(1971,2000)
    if(grepl("ff",mt)) it2 <- c(2071,2100) 
    if(grepl("nf",mt)) it2 <- c(2021,2050)
    label.title <- paste(attr(ceof[[2]],"longname")," change\n",
                         "ensemble ",FUN," (",
                         paste(it1,collapse="-")," to ",paste(it2,collapse="-"),")",sep="")
  } else if(grepl("mv",mt)) {
    if(is.null(colbar$pal)) {
      if(grepl("t2m|tas|temp",attr(ceof[[2]],"variable"))) {
        if(is.null(colbar$pal)) colbar$pal <- "t2m"
        colbar$breaks <- seq(-60,40,5)
      }
      if(grepl("pr",attr(ceof[[2]],"variable"))) {
        if(is.null(colbar$pal)) colbar$pal <- "precip"
        colbar$breaks <- seq(0,15,0.5)
      }
    }
    if(grepl("ff",mt)) it1 <- c(2071,2100) 
    if(grepl("nf",mt)) it1 <- c(2021,2050)
    if(grepl("pd",mt)) it1 <- c(1971,2000)
    it2 <- NULL
    label.title <- paste("Ensemble mean of ",attr(ceof[[2]],"longname")," (",paste(it1,collapse="-"),")",sep="")
  }
  Y1 <- map.commonEOF(ceof,is=is,im=im,ip=ip,it=it1,FUN=FUN,plot=FALSE)
  Y <- apply(Y1,2,mean,na.rm=TRUE)
  if(!is.null(it2)) {
    Y2 <- map.commonEOF(ceof,is=is,im=im,ip=ip,it=it2,FUN=FUN,plot=FALSE)
    dY <- apply(Y2,2,mean,na.rm=TRUE) - Y
    Y <- dY
  } 
  dim(Y) <- c(1,length(Y))
  Y <- as.field(Y,1,attr(Y1,"longitude"),attr(Y1,"latitude"),
                param=attr(Y1,"variable"),unit=attr(Y1,"unit"))
  if(is.null(colbar$breaks))  {
    if(grepl("cc",mt)) {
      colbar$breaks <- pretty(c(-max(abs(Y),na.rm=TRUE),max(abs(Y),na.rm=TRUE)),n=10)
    } else {
      colbar$breaks <- pretty(range(Y,na.rm=TRUE),n=10)
    }
  }
  map(Y,new=new,colbar=colbar,main=label.title)
  invisible(Y)
} 

scatterplot <- function(x,y,ix=NULL,xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,
                        main=NULL,legend=NULL,pal="cat",pch=21,cex=1.5,lwd=1.5,
                        new=FALSE,verbose=FALSE) {
  if(verbose) print("scatterplot")
  if(is.null(xlab)) xlab <- paste(attr(x,"variable")," (",attr(x,"unit"),")",sep="")
  if(is.null(ylab)) ylab <- paste(attr(y,"variable")," (",attr(y,"unit"),")",sep="")
  if(is.null(main)) main <- ""
  if(is.null(legend)) legend <- names(x)
  if(is.null(xlim)) xlim <- range(x,na.rm=TRUE) + c(-1,1)*diff(range(x,na.rm=TRUE))*0.1 
  if(is.null(ylim)) ylim <- range(y,na.rm=TRUE) + c(-1,1)*diff(range(y,na.rm=TRUE))*0.1
  if(!is.null(pal)) {
    col <- colscal(n=length(dtas),col=pal)
  } else {
    col <- rep("grey50",length(dtas))
  }
  bg <- col
  if(!is.null(ix)) col[im] <- "black"
  if(!is.null(ix)) {
    if(length(cex)==1) cex <- rep(cex,length(x))
    cex[ix] <- cex[1]*1.5
  }
  if(new) dev.new()
  plot(unlist(x),unlist(y),col=col,pch=pch,cex=cex,lwd=lwd,
       bg=bg,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  lines(xlim*1.5,rep(0,2),lwd=0.2)
  lines(rep(0,2),ylim*1.5,lwd=0.2)
  grid()
  legend("bottomleft",ncol=floor(length(x)/3),pch=pch,cex=1,col=col,pt.bg=bg,
         bg=adjustcolor("white",alpha=0.6),box.lwd=0.5,legend=legend)
}

select.colbar <- function(x,breaks=NULL,pal=NULL) {
  if(is.null(breaks)) {
    if(attr(x,"variable")=="correlation") {
      breaks <- seq(-1,1,0.1)
    } else if (attr(x,"unit") %in% c("%","pct")) {
      breaks <- seq(0,max(x,100,na.rm=TRUE),5)
    } else if(min(x,na.rm=TRUE)<0 & max(x,na.rm=TRUE)>0) {
      breaks <- pretty(c(-max(abs(x)*1.05,na.rm=TRUE),max(abs(x)*1.05,na.rm=TRUE)),n=8)
    } else if(min(x,na.rm=TRUE)>0 & max(x,na.rm=TRUE)>0) {
      breaks <- pretty(c(0,max(abs(x)*1.05,na.rm=TRUE)),n=8)
    } else if(min(x,na.rm=TRUE)<0 & max(x,na.rm=TRUE)<0) {
      breaks <- pretty(c(0,max(abs(x)*1.05,na.rm=TRUE)),n=8)
    }
  }
  if(is.null(pal)) {
    if(max(breaks)>0 & min(breaks)<0) {
      pal <- "burd"
    } else {
      pal <- "bu"
    }
  }
  return(list(pal=pal,breaks=breaks))
}
