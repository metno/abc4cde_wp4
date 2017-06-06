
# Get coordinates of a chosen prudence region
# Data about the regions are included in the folder back-end/data/PRUDENCE 
getPrudenceCoords <- function(prudence,region,destfile="coords.txt"){
  i <- region
  if(is.character(region)) i <- which(prudence[,2]==region)
  lat <- as.numeric(prudence[i,6:7])
  lon <- as.numeric(prudence[i,4:5])
  coords1 <- expand.grid(lon[1],lat)
  coords2 <- expand.grid(lon[2],lat)
  coords <- rbind(coords1,coords2[c(2,1),])
  write(t(coords),file=destfile,ncolumns = 2)
}

# Transform the coordinates of a region from Cartesian to polar. 
# There is most likely a version of this in esd. 
getPolCoords <- function(region,shape="",destfile="coords.txt"){
  if(is.character(region)) region <- which(as.character(shape$LAB)==region)
  pol.coords <- coordinates(shape@polygons[[region]]@Polygons[[1]])
  write(t(pol.coords),file=destfile,ncolumns = 2)
}

#apply mask to a zoo object by setting values outside the mask to NA
mask.zoo <- function(zoo.object,mask){
  mask <- flip(mask,direction='y')
  zoo.object[,which(is.na(getValues(mask)))] <- NA
  return(zoo.object)
}

#Get grid boxes belonging to a SREX region and calculate some basic statistics for it.
get.srex.region <- function(destfile,region=NULL,print.srex=FALSE,verbose=FALSE) {
  if(verbose) print("get.srex.region")
  home <- system("echo $HOME",intern=TRUE)
  shape <-  get.shapefile("referenceRegions.shp")
  X <- retrieve(destfile,lon=NULL,lat=NULL,verbose=verbose)
  srex <- list()
  if(is.null(region)){
    for (i in 1:length(levels(shape$LAB))){
      polygon <- shape[i,]
      mask <- gen.mask.srex(destfile=destfile, mask=polygon, 
                            ind=FALSE, inverse=FALSE, mask.values=1)
      if(verbose){
        if(i==1){
          plot(shape)
        }
        plot.mask <- mask
        extent(plot.mask) <- c(-180,180,-90,90)
        projection(plot.mask) <- projection(shape)
        plot(plot.mask,col=rainbow(100, alpha=0.35)[sample(1:100,1)],legend=F,add=T)
      }
      name <- levels(shape$NAME)[i]
      X.region <- mask.zoo(X,mask)
      srex[[name]]$name <- name
      srex[[name]]$label <- levels(shape$LAB)[i]
      srex[[name]]$area.mean <- aggregate.area(X.region,FUN="mean",na.rm=T)
      srex[[name]]$area.sd <- aggregate.area(X.region,FUN="sd",na.rm=T)
    }  
  } else {
    polygon <- shape[levels(shape$LAB)==region,]
    mask <- gen.mask.srex(destfile=destfile, mask=polygon, ind=FALSE, 
                          inverse=FALSE, mask.values=1)
    if(verbose) {
      plot(shape)
      plot.mask <- mask
      extent(plot.mask) <- c(-180,180,-90,90)
      projection(plot.mask) <- projection(shape)
      plot(plot.mask,col=rainbow(100, alpha=0.35)[sample(1:100,1)],legend=F,add=T)
    }
    name <- levels(shape$NAME)[i]
    X.region <- mask.zoo(X,mask)
    srex[[name]]$name <- name
    srex[[name]]$label <- levels(shape$LAB)[i]
    srex[[name]]$area.mean <- aggregate.area(X.region,FUN="mean",na.rm=T)
    srex[[name]]$area.sd <- aggregate.area(X.region,FUN="sd",na.rm=T) 
  }
  
  if(print.srex) {
    print("Region names in alphabetical order and the corresponding label to be used when selecting the region:")
    print(data.frame(NAME=gsub("\\[[^\\]]*\\]", "", levels(shape$NAME), perl=TRUE),
                     LABEL=levels(shape$LAB)))
    return()
  }
  return(srex)
}

# Create a raster mask for the selected SREX sub-region from the CMIP5 netcdf file.
gen.mask.srex <- function(destfile, mask.polygon=NULL, ind=FALSE, inverse=FALSE, mask.values=1){
  print(destfile)
  r <- raster(destfile)
  r <- setValues(r,NA)
  extent.r <- extent(r)
  if(extent.r[2]==360) extent(r) <- c(-180,180,-90,90)
  indices <- extract(r,mask.polygon,cellnumbers=TRUE)[[1]][,1]
  if(extent(mask.polygon)[2]>180){
    extent(r) <- c(180,540,-90,90)
    indices <- sort(c(indices,extract(r,mask.polygon,cellnumbers=TRUE)[[1]][,1]))
  }
  if(inverse){
    tmp <- seq(1,length(getValues(r)))
    indices <- tmp[which(is.na(match(tmp,indices)))]
  }
  mask.raster <- r
  extent(mask.raster) <- c(0,360,-90,90)
  mask.raster[indices] <- mask.values
  if(ind) return(indices)
  return(mask.raster)
}
