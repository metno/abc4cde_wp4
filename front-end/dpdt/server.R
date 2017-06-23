
library(shiny)
library(DECM)
library(DT)

## load statistics
stats <- NULL
data("statistics.cmip.era.tas.1981-2010")
stats$tas$present <- store
data("statistics.cmip.tas.2021-2050")
stats$tas$nf <- store
data("statistics.cmip.tas.2071-2100")
stats$tas$ff <- store
data("statistics.cmip.era.pr.1981-2010")
stats$pr$present <- store
data("statistics.cmip.pr.2021-2050")
stats$pr$nf <- store
data("statistics.cmip.pr.2071-2100")
stats$pr$ff <- store

dT <- list()
dPr <- list()
gcms <- names(stats$tas$ff)
for (gcm in gcms) {
  dT[["nf"]][[gcm]] <- stats$tas$nf[[gcm]]$mean - stats$tas$present[[gcm]]$mean
  dT[["ff"]][[gcm]] <- stats$tas$ff[[gcm]]$mean - stats$tas$present[[gcm]]$mean
  dPr[["nf"]][[gcm]] <- stats$pr$nf[[gcm]]$mean - stats$pr$present[[gcm]]$mean
  dPr[["ff"]][[gcm]] <- stats$pr$ff[[gcm]]$mean - stats$pr$present[[gcm]]$mean
}

# Define a server for the Shiny app
shinyServer(function(input, output) {
  
  output$dtdpr <- renderPlot({
    season <- switch(tolower(as.character(input$season)),
                     'annual mean'='ann','winter'='djf','spring'='mam',
                     'summer'='jja','autumn'='son')
    period <- switch(tolower(as.character(input$period)),
                     "far future (2071-2100)"='ff',
                     "near future (2021-2050)"='nf')
    dtas <- unlist(lapply(dT[[period]], function(x) x[[season]]))
    dpr <- unlist(lapply(dPr[[period]], function(x) x[[season]]))
    scatterplot(dtas,dpr,ix=NULL,xlim=c(-3,3),ylim=c(-0.2,0.2)/(60*60*24),
                xlab="Temperature (deg C)",ylab="Precipitation (kg m-1 s-2)",
                main=paste("Climate change assuming RCP4.5\npresent day to",input$period),
                show.legend=FALSE,
                legend=seq(length(dtas)),pal=NULL,#pal="cat",pch=21,
                pch=as.character(seq(length(dtas))),cex=1.5,lwd=1.5,new=FALSE)
  }, width=600, height=600)
  
})