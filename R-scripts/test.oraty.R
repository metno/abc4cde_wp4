## test.ratyo.R
## Test Olle's back-end functions
library(DECM)

file.meta <- find.file("metaextract.rda")
file.meta <- file.meta[grep("abc4cde",file.meta)]

file.shape <- find.file("*.shp")
file.shape <- file.shape[grep("abc4cde",file.shape)]
X <- get.shapefile(file.shape,with.path=TRUE)

## Compare a) cdo.mean b) reading the data into R and computing the mean value c) calculating the mean from the commmonEOFs
## a & b)
variable <- "pr"
i <- 1
it <- c(2071,2100)#c(1981,2010)
X <- getGCMs(select=i,varid=variable)
print(paste("unit in file:",X[[1]]$var$pr$units))
gcm.file <- X[[1]]$filename
m.a <- c(cdo.mean(gcm.file,it),cdo.mean(gcm.file,it,seasonal=T))
Y <- retrieve(gcm.file,it=it)
print(paste("unit after esd retrieve:",attr(Y,"unit")))
m.b <- c(mean(aggregate.area(Y,FUN="mean")),
         mean(aggregate.area(subset(Y,it="djf"),FUN="mean")),mean(aggregate.area(subset(Y,it="mam"),FUN="mean")),
         mean(aggregate.area(subset(Y,it="jja"),FUN="mean")),mean(aggregate.area(subset(Y,it="son"),FUN="mean")))
m.a*(60*60*24)

data("ceof.gcm.pr.annual")
Z <- map.commonEOF(ceof,im=i,it=it)
m.c <- mean(aggregate.area(Z,FUN="mean"))

# missing script getMonthlyERA.py
#getERA("tas",start=1979,end=2016,griddes="cmip_1.25deg_to_2.5deg.txt",destfile="test.nc")
#python.getEra(1979,2016,"tas",0,type,stream,outfile)

getCFSR()

verbose <- TRUE
X <- getGCMs(select=1:9,varid="tas",verbose=verbose)
Y <- getGCMs(select=1:9,varid="pr",verbose=verbose)
X.all$mean <- lapply(X,function(x) cdo.mean(x,period=c(1850,2100)))
X.all$spatSd <- lapply(X,function(x) cdo.spatSd(x,period=c(1850,2100)))
X.present$mean <- lapply(X,function(x) cdo.mean(x,period=c(1971,2000)))
X.present$spatSd <- lapply(X,function(x) cdo.spatSd(x,period=c(1971,2000)))
#X.present$timeSd <- lapply(X,function(x) cdo.timeSd(x,period=c(1971,2000)))
X.near$mean <- lapply(X,function(x) cdo.mean(x,period=c(2021,2050)))
X.near$spatSd <- lapply(X,function(x) cdo.spatSd(x,period=c(2021,2050)))
X.far$mean <- lapply(X,function(x) cdo.mean(x,period=c(2071,2100)))
X.far$spatSd <- lapply(X,function(x) cdo.spatSd(x,period=c(2071,2100)))

Y.all$mean <- lapply(Y,function(x) cdo.mean(x,period=c(1850,2100)))
Y.all$spatSd <- lapply(Y,function(x) cdo.spatSd(x,period=c(1850,2100)))
Y.present$mean <- lapply(Y,function(x) cdo.mean(x,period=c(1971,2000)))
Y.present$spatSd <- lapply(Y,function(x) cdo.spatSd(x,period=c(1971,2000)))
#Y.present$timeSd <- lapply(Y,function(x) cdo.timeSd(x,period=c(1971,2000)))
Y.near$mean <- lapply(Y,function(x) cdo.mean(x,period=c(2021,2050)))
Y.near$spatSd <- lapply(Y,function(x) cdo.spatSd(x,period=c(2021,2050)))
Y.far$mean <- lapply(Y,function(x) cdo.mean(x,period=c(2071,2100)))
Y.far$spatSd <- lapply(Y,function(x) cdo.spatSd(x,period=c(2071,2100)))
