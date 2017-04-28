# plottingtools.R

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
      pal <- "t2m"
    } else {
      pal <- "bu"
    }
  }
  return(list(pal=pal,breaks=breaks))
}
