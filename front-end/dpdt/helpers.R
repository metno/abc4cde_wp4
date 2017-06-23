## helpers.R
## Help functions for the shiny app "dpdt"

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