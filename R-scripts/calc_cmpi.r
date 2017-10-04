library(DECM)
source("/home/ubuntu/abc4cde_wp4/back-end/R/") #Did not manage to install the latesr version of DECM!

#Calculate weights for the weighted average of rms, assuming monthly time step.
calc.mon.weights <- function(lon,lat){
  
  weights <- array(NA,dim=c(12,length(lon),length(lat)))
  time.weights <- c(31,28,31,30,31,30,31,31,30,31,30,31)/365
  lat.weights <- rep(cos(pi*lat/180),length(lon))
  dim(lat.weights) <- c(length(lon),length(lat))
  lon.weights <- rep(1,length(lon)*length(lat))
  dim(lon.weights) <- c(length(lon),length(lat))
  
  for(i in 1:length(lat)){
    weights[,,i] <- time.weights
  }
  
  for(i in 1:12){
    weights[i,,] <- weights[i,,]*lat.weights
  }
  
  return(weights)
}

period <- c(1981,2010)
shape <-  get.shapefile("/home/ubuntu/abc4cde_wp4/back-end/data/SREX_regions/referenceRegions.shp",with.path = T)
 getGCMs(select=1:4,varid="pr")
 getGCMs(select=1:4,varid="tas")

vars <- c("tas","pr")

#Pre-process ERA-interim
for(var in vars){
  era.file <- paste("/home/ubuntu/Data/Data/era-interim_monthly_1979-2017_",var,".2.5deg.nc",sep="")
  era.mulc <- paste("/home/ubuntu/Data/Data/era.mulc",var,"nc",sep=".")
  era.mon.file <- paste("/home/ubuntu/Data/Data/era.monmean",var,"nc",sep=".")
  if(var=="pr"){
    cdo.command("mulc",1000,era.file,era.mulc)
  }else{
    era.mulc <- era.file
  }
  cdo.command(c("-ymonmean","-selyear"),c("",paste(period,collapse="/")),
              infile=era.mulc,outfile=era.mon.file)
}

#Calculate weights only once
r <- raster(era.mon.file)
lon <- xFromCol(r)
lat <- yFromRow(r)
weights <- calc.mon.weights(lon,lat)

#For testing
vars <- c("tas","pr")
rms <- e <- list()
models <- 1:4
for(var in vars){
  era.mon.file <- paste("/home/ubuntu/Data/Data/era.monmean",var,"nc",sep=".")
  era <- coredata(retrieve(era.mon.file))
  dim(era) <- c(12,length(longitude(era)),length(latitude(era)))
  for(i in models){
    gcm.file <- paste(paste("/home/ubuntu/Data/Data/GCM",i,sep=""),var,"nc",sep=".")
    gcm.mon.file <- "/home/ubuntu/Data/Data/gcm.monmean.nc"
    cdo.command(c("-ymonmean","-selyear"),c("",paste(period,collapse="/")),
                infile=gcm.file,outfile=gcm.mon.file)
    gcm <- coredata(retrieve(gcm.mon.file))
    dim(gcm) <- c(12,length(longitude(gcm)),length(latitude(gcm)))
    rms[[var]][[i]] <- sqrt(sum(weights*(gcm-era)^2)/sum(weights))
  }
  
  rms[[var]][["median"]] <- median(unlist(rms[[var]]))
  for(i in 1:4){
    e[[var]][[i]] <- (rms[[var]][[i]]-rms[[var]][["median"]])/rms[[var]][["median"]]
  }
}
