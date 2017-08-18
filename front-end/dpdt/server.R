
library(shiny)
library(DECM)
library(DT)
source("helpers.R")

## Load statistics calculated with script 'calculate_statistics.R'
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

## Function 'regions' is defined in helpers.R
srex <- regions("srex")

## Load geographical data for map
data("geoborders",envir=environment())
  
## Define a server for the Shiny app
shinyServer(function(input, output) {
  
  output$dtdpr <- renderPlot({
    season <- switch(tolower(as.character(input$season)),
                     'annual mean'='ann',
                     'winter'=c('dec','jan','feb'),
                     'spring'=c('mar','apr','may'),
                     'summer'=c('jun','jul','aug'),
                     'autumn'=c('sep','oct','nov'))
                     #'annual mean'='ann','winter'='djf','spring'='mam',
                     #'summer'='jja','autumn'='son')
    period <- switch(tolower(as.character(input$period)),
                     "far future (2071-2100)"='ff',
                     "near future (2021-2050)"='nf')
    gcms <- names(stats$tas$ff)
    if(tolower(input$region)=="global") {
      coord <- list(lon=c(-180,180),lat=c(-90,90))
      dtas <- sapply(gcms, function(gcm) mean(sapply(season, function(s)
                       stats$tas[[period]][[gcm]][["mean"]][[s]])) - 
                       mean(sapply(season, function(s)
                       stats$tas$present[[gcm]][["mean"]][[s]])))
      dpr <- sapply(gcms, function(gcm) mean(sapply(season, function(s)
                       stats$pr[[period]][[gcm]][["mean"]][[s]])) - 
                       mean(sapply(season, function(s)
                       stats$pr$present[[gcm]][["mean"]][[s]])))
      #dtas <- sapply(gcms, function(gcm) stats$tas[[period]][[gcm]][["mean"]][[season]] - 
      #                 stats$tas$present[[gcm]][["mean"]][[season]]) 
      #dpr <- sapply(gcms, function(gcm) stats$pr[[period]][[gcm]][["mean"]][[season]] - 
      #                stats$pr$present[[gcm]][["mean"]][[season]])
    } else {
      i.srex <- which(srex$name==input$region)
      region <- srex$label[i.srex]
      dtas <- sapply(gcms, function(gcm) mean(sapply(season, function(s)
                       stats$tas[[period]][[gcm]][[region]][["mean"]][[s]])) - 
                       mean(sapply(season, function(s) 
                       stats$tas$present[[gcm]][[region]][["mean"]][[s]])))
      dpr <- sapply(gcms, function(gcm) mean(sapply(season, function(s)
                       stats$pr[[period]][[gcm]][[region]][["mean"]][[s]])) - 
                       mean(sapply(season, function(s) 
                       stats$pr$present[[gcm]][[region]][["mean"]][[s]])))
      #dtas <- sapply(gcms, function(gcm) stats$tas[[period]][[gcm]][[region]][["mean"]][[season]] - 
      #                stats$tas$present[[gcm]][[region]][["mean"]][[season]])
      #dpr <- sapply(gcms, function(gcm) stats$pr[[period]][[gcm]][[region]][["mean"]][[season]] - 
      #                                  stats$pr$present[[gcm]][[region]][["mean"]][[season]])
    }
    im <- as.numeric(gsub(":.*","",input$gcms))
    scatterplot(dtas,dpr*(60*60*24),ix=NULL,xlim=input$tlim,ylim=input$plim,
                xlab="Temperature change (deg C)",ylab="Precipitation change (mm/day)",
                main=paste("Climate change assuming RCP4.5\npresent day (1981-2010) to",input$period),
                show.legend=FALSE,im=im,
                legend=seq(length(dtas)),pal=NULL,#pal="cat",pch=21,
                pch=as.character(seq(length(dtas))),cex=1.5,lwd=1.5,new=FALSE)
  }, width=function(){450}, height=function(){450})
  
  output$map <- renderPlot({
    if(tolower(input$region)=="global") {
      region <- list(lon=c(-180,-180,180,180,-180),lat=c(-90,90,90,-90,-90))
    } else {
      i.srex <- which(srex$name==input$region)
      region <- list(lon=srex$coord[[i.srex]][1,],
                     lat=srex$coord[[i.srex]][2,])
    }
    par(mgp=c(1,0.5,0),mar=c(0.2,0.2,0.2,0.2))
    plot(geoborders$x,geoborders$y,col="grey30",type="l",lwd=0.5,
         xlim=c(-180,180),ylim=c(-90,90),
         xlab="Longitude",ylab="Latitude",xaxt="n",yaxt="n")
    #lines(attr(geoborders,'borders')$x,attr(geoborders,'borders')$y,col="grey30")
    par(xaxt="s",yaxt="s",las=1,col.axis='grey',col.lab='grey20',
        cex.lab=0.7,cex.axis=0.7)
    axis(3,at=pretty(par("xaxp")[1:2],n=5),col='grey50')
    axis(2,at=pretty(par("yaxp")[1:2],n=5),col='grey50')
    grid()
    lines(region$lon,region$lat,col="blue",lwd=1.5,lty=1)
  }, width=function(){200},height=function(){200}*0.6)#width=250, height=175)
  
})