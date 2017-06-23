## helpers.R
## Help functions for the shiny app "modelexplorer"

map.changes <- function(ceof,period="ff",FUN="mean",new=FALSE) {
  if(new) dev.new()
  par(fig=c(0,0.5,0,0.33), new=FALSE)  
  map.ensemble(ceof,type="mvpd",new=FALSE,FUN=FUN)
  par(fig=c(0.5,1,0,0.33), new=FALSE)  
  map.ensemble(ceof,type=paste("mv",period,sep=""),new=FALSE,FUN=FUN)
  par(fig=c(0,1,0.33,1), new=FALSE)  
  map.ensemble(ceof,type=paste("cc",period,sep=""),new=FALSE,FUN=FUN)
}

fn.switch <- function(fn="mean") {
  if(!is.null(fn)) {
    switch(fn,"Mean"="mean",
           "95th percentile"="q95",
           "5th percentile"="q5")
  } else {
    "mean"
  }
}

type.switch <- function(type="Climate change far future") {
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