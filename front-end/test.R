# test.R

## Set paths so that R can find data and R-scripts
path.data <- "back-end/data"
path.r <- "back-end/R"

source(file.path(path.r,"cds.R"))
source(file.path(path.r,"eqc.R"))

## Compare fields obtained from common EOF to original data
varid <- "pr"
it <- "annual"
is <- list(lon=c(-180,180),lat=c(0,90))
label.is <- paste(paste(sapply(is$lon, function(x) if(x<0) paste(abs(x),"W",sep="") else paste(x,"E",sep="")),collapse="-"), 
            paste(sapply(is$lat, function(x) if(x<0) paste(abs(x),"S",sep="") else paste(x,"N",sep="")),collapse="-"),sep=".")

## Load and prepare GCM commonEOF
fname <- file.path(path.data,paste("ceof.gcm",varid,it,sep="."))
if(!is.null(is)) fname <- paste(fname,label.is,sep=".")
fname <- paste(fname,"rda",sep=".")
if(!file.exists(fname)) {
  destfile <- paste(path.data,"/GCM",1:9,".",varid,".nc",sep="")
  ceof <- commonEOFS.gcm(varid=varid,it=it,is=is,destfile=destfile,destfile.ceof=fname,verbose=TRUE)
} else {
  load(fname)
}
Y <- as.field.commonEOFS(ceof)

## Load and prepare RCM commonEOF
fname <- file.path(path.data,paste("ceof.rcm",varid,it,"rda",sep="."))
if(!file.exists(fname)) {
  destfile <- paste(path.data,"/CM",1:9,".",varid,".nc",sep="")
  ceof <- commonEOFS.rcm(varid=varid,it=it,destfile=destfile,destfile.ceof=fname,verbose=TRUE)
} else {
  load(fname)
}
Yb <- as.field.commonEOFS(ceof)

i <- 5
f <- getGCMs(select=i,varid=varid,destfile=file.path(path.data,paste("GCM",i,".",varid,".nc",sep="")))
x <- retrieve(f[[1]]$filename)
if(!is.null(is)) x <- subset(x,is=is)
if(it=="annual") x <- annual(x) else x <- subset(x,it=it)
y <- Y[[i]]
z <- compare.fields(x,y,lplot=FALSE,verbose=TRUE)
fname <- paste("eval.commoneofs.gcm",i,".",varid,".",it,sep="")
if(!is.null(is)) fname <- paste(fname,label.is,sep=".")
compare.fields(z,lplot=TRUE,type=c("correlation","rmsd"),
  filename=paste(fname,"corr.rmsd.pdf",sep="."))
compare.fields(z,lplot=TRUE,type=c("meandiff","trenddiff"),
  filename=paste(fname,"mean.trend.pdf",sep="."))

f <- getRCMs(select=i,varid=varid,destfile=file.path(path.data,paste("CM",i,".",varid,".nc",sep="")))
xb <- retrieve(f[[1]]$filename)
if(it=="annual") xb <- annual(xb) else xb <- subset(xb,it=it)
yb <- Yb[[i]]
zb <- compare.fields(xb,yb,lplot=FALSE,verbose=TRUE)
fname <- paste("eval.commoneofs.rcm",i,".",varid,".",it,sep="")
compare.fields(zb,lplot=TRUE,type=c("correlation","rmsd"),
  filename=paste(fname,"corr.rmsd.pdf",sep="."))
compare.fields(zb,lplot=TRUE,type=c("meandiff","trenddiff"),
  filename=paste(fname,"mean.trend.pdf",sep="."))


