# test.R

## Set paths so that R can find data and R-scripts
path.data <- "back-end/data"
path.r <- "back-end/R"

source(file.path(path.r,"cds.R"))
source(file.path(path.r,"eqc.R"))

varid <- "tas"

## Compare GCM field obtained from common EOF to original data
load(file.path(path.data,paste("ceof.gcm",varid,"annual.rda",sep=".")))
Y <- as.field.commonEOFS(ceof)

i <- 2
fn <- getGCMs(select=i,varid="tas")
x <- annual(retrieve(fn[[1]]$filename))
y <- Y[[i]]
z <- compare.fields(x,y,lplot=TRUE)

## Compare RCM field obtained from common EOF to original data
load(file.path(path.data,paste("ceof.rcm",varid,"annual.rda",sep=".")))
Yb <- as.field.commonEOFS(ceof)

i <- 2
fn <- getRCMs(select=i,varid="tas")
xb <- annual(retrieve(fn[[1]]$filename))
yb <- Yb[[i]]
zb <- compare.fields(xb,yb,lplot=FALSE)


