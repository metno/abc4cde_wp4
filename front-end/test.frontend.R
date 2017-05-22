## test.frontend.R
## test some ideas for the front-end
## author: Kajsa Parding
## created: 11 May 2017
## last updated: 11 May 2017

path.data <- "/home/kajsamp/git/abc4cde_wp4/back-end/data"
path.R <- "/home/kajsamp/git/abc4cde_wp4/back-end/R"
source(file.path(path.R,"cds.R"))
load(file.path(path.data,"ceof.gcm.tas.annual.rda"))
ceof.tas.gcms <- ceof
load(file.path(path.data,"ceof.gcm.pr.annual.rda"))
ceof.pr.gcms <- ceof
tas.gcms <- as.field.commonEOFS(ceof.tas.gcms)
pr.gcms <- as.field.commonEOFS(ceof.pr.gcms)

dtdprplot <- function(tas,pr=NULL,im=c(2),pal="cat",lplot=TRUE,new=TRUE) {

  dt <- function(x,it1=c(1981,2010),it2=c(2071,2100),is=NULL,FUN="mean") {
    if(!is.null(is)) x <- subset(x,is=is)
    x1 <- mean(aggregate.area(subset(x,it=it1),FUN=FUN))
    x2 <- mean(aggregate.area(subset(x,it=it2),FUN=FUN))
    return(x2-x1)
  }

  if(is.null(pr) & length(tas)==2) {
    dtas <- tas[[1]]
    dpr <- tas[[2]]
  } else {
    if(is.null(dim(tas[[1]]))) {
      dtas <- tas
    } else {
      dtas <- lapply(tas,dt)
    }
    if(is.null(dim(tas[[1]]))) {
      dpr <- pr
    } else {
      dpr <- lapply(pr,dt)  
    }
  }
  
  if(lplot) {
    xlim <- c(-3,3)
    ylim <- c(-0.2,0.2)
    if(!is.null(pal)) {
      col <- colscal(n=length(dtas),col=pal)
    } else {
      col <- rep("grey50",length(dtas))
    }
    col[im] <- "black"
    pch <- 1:length(dtas)
    cex <- rep(1,length(dtas))
    cex[im] <- 2.5
    lwd <- rep(1.5,length(dtas))
    lwd[im] <- 2.5
    ## Add rcp label  
    if(new) dev.new()
    plot(unlist(dtas),unlist(dpr),col=col,pch=pch,cex=cex,lwd=lwd,
         xlim=xlim,ylim=ylim,
         xlab="Temperature change (degC)",
         ylab="Precipitation change (mm/day)",
         main="Global mean climate change\n(1981-2010 to 2071-2100)")
    lines(xlim*1.5,rep(0,2),lwd=0.2)
    lines(rep(0,2),ylim*1.5,lwd=0.2)
    grid()
    legend("bottomleft",ncol=2,pch=pch,cex=0.8,col=col,
           bg=adjustcolor("white",alpha=0.6),box.lwd=0.5,
           legend=paste(attr(tas,"model_id")$gcm,
                        attr(tas,"model_id")$gcm_rip,sep="."))
  }

  dtdpr <- list(dtas=dtas,dpr=dpr)
  attr(dtdpr,"model_id") <- attr(tas,"model_id")
  return(dtdpr)
}

dtdpr <- dtdprplot(tas.gcms,pr.gcms)
dtdprplot(dtdpr,pal="rainbow")

