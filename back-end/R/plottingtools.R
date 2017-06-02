
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
    if(is.null(colbar$pal)) colbar$pal <- "burd"
    it1 <- c(1971,2000)
    if(grepl("ff",mt)) it2 <- c(2071,2100) 
    if(grepl("nf",mt)) it2 <- c(2021,2050)
    label.title <- paste(attr(ceof[[2]],"longname")," change\n",
                         "ensemble ",FUN," (",
                         paste(it1,collapse="-")," to ",paste(it2,collapse="-"),")",sep="")
  } else if(grepl("mv",mt)) {
    if(is.null(colbar$pal)) {
      if(grepl("t2m|tas|temp",attr(ceof[[2]],"variable"))) colbar$pal <- "t2m"
      if(grepl("pr",attr(ceof[[2]],"variable"))) colbar$pal <- "precip"
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
