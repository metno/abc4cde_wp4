## Kajsa Parding, 2017-05-15
## R-shiny app that presents GCM and RCM data.

library(shiny)
library(DECM)
library(DT)
source("helpers.R")
#if ('RgoogleMaps' %in% installed.packages()) install.packages('RgoogleMaps')

## Preparations
## source scripts
## load metadata
data("metaextract")
M <- data.frame(list(project_id=meta$project_id, experiment_id=meta$experiment_id, gcm=meta$gcm,
                     rip=meta$gcm_rip, rcm=meta$rcm, var=meta$var, unit=meta$unit, resolution=paste(meta$resolution,"deg"),
                     domain=paste(gsub(","," - ",meta$lon),"E"," / ",paste(gsub(","," - ",meta$lat)),"N",sep=""), 
                     years=gsub(",","-",gsub("-[0-9]{2}","",meta$dates)), url=meta$url))
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

## load and expand commonEOFs
ceof.all <- NULL
data("ceof.gcm.tas.annual")
ceof.all$tas$CMIP5 <- ceof
data("ceof.gcm.pr.annual")
ceof.all$pr$CMIP5 <- ceof
data("ceof.rcm.tas.annual")
ceof.all$tas$CORDEX <- ceof
data("ceof.rcm.pr.annual")
ceof.all$pr$CORDEX <- ceof
ceof <- ceof.all
rm("ceof.all"); gc(reset=TRUE)
selectrowindex <- 1

dT <- list()
dPr <- list()
gcms <- names(stats$tas$ff)
for (gcm in gcms) {
  dT[["nf"]][[gcm]] <- stats$tas$nf[[gcm]]$mean - stats$tas$present[[gcm]]$mean
  dT[["ff"]][[gcm]] <- stats$tas$ff[[gcm]]$mean - stats$tas$present[[gcm]]$mean
  dPr[["nf"]][[gcm]] <- stats$pr$nf[[gcm]]$mean - stats$pr$present[[gcm]]$mean
  dPr[["ff"]][[gcm]] <- stats$pr$ff[[gcm]]$mean - stats$pr$present[[gcm]]$mean
}

select.ceof <- function(table_rows_selected=1,varid="Temperature") {
  if (length(table_rows_selected)>0) {
    selectedrowindex <- table_rows_selected[length(table_rows_selected)]
    selectedrowindex <- as.numeric(selectedrowindex)
  } else {
    selectedrowindex <- 1
  }
  selectedrow <- (M[selectedrowindex,])
  ceof.sel <- NULL
  if (grepl("temp",tolower(varid))) {
    ceof.sel <- ceof$tas
  } else {
    ceof.sel <- ceof$pr
  }
  ceof.sel <- ceof.sel[[selectedrow$project_id]]
  im <- which(attr(ceof,"model_id")$rcm==selectedrow$rcm & 
              attr(ceof,"model_id")$gcm==selectedrow$gcm & 
              attr(ceof,"model_id")$gcm_rip==selectedrow$rip)
  attr(ceof.sel,"im") <- im
  return(ceof.sel)
}

# Define a server for the Shiny app
shinyServer(function(input, output) {
  
  # Filter data based on selections
  output$table <- DT::renderDataTable({
    datatable(M,filter='top',selection = 'single',options = list(), style="bootstrap")
  })
  
  output$text <- renderText({
    txt <- "Filter your search and click on one climate model"
  })
  
  output$glossary <- DT::renderDataTable({
    datatable(rbind(c('tas','Near-surface air temperature'),
                    c('pr','Precipitation'),
                    c('GCM','Global Climate Model or General Circulation Model'),
                    c('rip','Realisation, initialisation method, and physics version of a GCM'),
                    c('RCM','Regional Climate Model')),
              colnames = c('',''),options = list(pageLength = 10), 
              caption='Glossary and variable names',style="bootstrap")
  })
  
  output$selectedrow <- DT::renderDataTable({
    if (length(input$table_rows_selected)>0) {
      selectedrowindex <- input$table_rows_selected[length(input$table_rows_selected)]
      selectedrowindex <- as.numeric(selectedrowindex)
    } else {
      selectedrowindex <- 1
    }
    selectedrow <- (M[selectedrowindex,])
    selectedrow
  })
  
  selected.ceof <- reactive({
    select.ceof(input$table_rows_selected,input$varid)  
  })
  
  output$dtdpr <- renderPlot({
    #ceof.tas <- select.ceof(input$table_rows_selected,"temp")
    #ceof.pr <- select.ceof(input$table_rows_selected,"precip")
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

  output$map <- renderPlot({
    ceof <- selected.ceof()
    fn <- 
    map.ensemble(ceof,type=type.switch(input$type),new=FALSE,FUN=fn.switch(input$fn))
  },height=500, width=750)
  
  output$rawdata <- DT::renderDataTable({
    #browser()
    # if (length(input$table_rows_selected)>0) {
    #   selectedrowindex <- input$table_rows_selected[length(input$table_rows_selected)]
    #   selectedrowindex <- as.numeric(selectedrowindex)
    # } else 
    #   selectedrowindex <- 1
    # selectedrow <- (select.station()[selectedrowindex,])
    # selectedrow
    # sel.sta <- station(select.station()[selectedrowindex,])  
    # # browser()
    # DT::datatable(data.frame(Date=index(sel.sta),Value=coredata(sel.sta)),filter='top',options = list(paging = TRUE))
  })
  
  
  # output$plots <- renderPlot({
  #   
  #   variable <- table[selectedrowindex,1]
  #   #write your plot function
  # })
  
  # downloadHandler() takes two arguments, both functions.
  # The content function is passed a filename as an argument, and
  #   it should write out data to that filename.
  
  #output$plot <- renderPlot({
  # library(esd)
  # data(Oslo)
  # # y <- station(datasetInput())
  # plot(Oslo,new=FALSE)
  #})
  
  output$downloadMeta <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    # filename = function() {
    #   library(esd)
    #   if (length(input$table_rows_selected)>0) {
    #     selectedrowindex <- input$table_rows_selected[length(input$table_rows_selected)]
    #     selectedrowindex <- as.numeric(selectedrowindex)
    #   } else
    #     selectedrowindex <- 1
    #   #browser()
    #   param <- ele2param(ele = select.station()[selectedrowindex,7],src=select.station()[selectedrowindex,10])
    #   loc <- select.station()[selectedrowindex,2]
    #   src <-  select.station()[selectedrowindex,10]
    #   cntr <- select.station()[selectedrowindex,3]
    #   paste(gsub(' ','-',param[2]),'_',loc,'-',cntr,'_',src,'meta.csv', sep = "")
    # },
    # 
    # # This function should write data to a file given to it by
    # # the argument 'file'.
    # content = function(file) {
    #   library(esd)
    #   if (length(input$table_rows_selected)>0) {
    #     selectedrowindex <- input$table_rows_selected[length(input$table_rows_selected)]
    #     selectedrowindex <- as.numeric(selectedrowindex)
    #   } else
    #     selectedrowindex <- 1
    #   # Write to a file specified by the 'file' argument
    #   write.csv(as.data.frame(select.station()[selectedrowindex,]), file, row.names = FALSE)
    # }
  )
  
  output$downloadData <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    # filename = function() {
    #   library(esd)
    #   if (length(input$table_rows_selected)>0) {
    #     selectedrowindex <- input$table_rows_selected[length(input$table_rows_selected)]
    #     selectedrowindex <- as.numeric(selectedrowindex)
    #   } else
    #     selectedrowindex <- 1
    #   param <- ele2param(ele = select.station()[selectedrowindex,7],src=select.station()[selectedrowindex,10])
    #   loc <- select.station()[selectedrowindex,2]
    #   src <-  select.station()[selectedrowindex,10]
    #   cntr <- select.station()[selectedrowindex,3]
    #   paste(gsub(' ','-',param[2]),'_',loc,'-',cntr,'_',src,'data.csv', sep = "")
    # },
    # 
    # # This function should write data to a file given to it by
    # # the argument 'file'.
    # content = function(file) {
    #   # Write to a file specified by the 'file' argument
    #   data <- esd::station(select.station()[1,])
    #   #save(data,file=file)
    #   write.table(data.frame(Date=index(data),Value=coredata(data)), file=file, row.names = FALSE)
    # }
  )
  
})


## Estimate the probabilities for trend in observation is within the population trends based on of downscaled results
## zoo objects are slow so extract the core data
# trendscore <- function(x) {
#   it.X <- year(x) 
#   X <- coredata(x)
#   it.y <- year(attr(x,'station'))
#   y <- coredata(attr(x,'station'))
#   X <- X[is.element(it.X,it.y),]
#   ty <- trend.coef(y)
#   tX <- apply(X,2,FUN='trend.coef')
#   score <- pnorm(ty,mean=mean(tX),sd=sd(tX))
#   return(c(score,lon(x),lat(x)))
# }
# 
# ## Estimate the probabilities for observed values are within the 90% conf. int. of population of downscaled results
# ## zoo objects are slow so extract the core data
# varscore <- function(x) {
#   it.X <- year(x) 
#   X <- coredata(x)
#   it.y <- year(attr(x,'station'))
#   y <- coredata(attr(x,'station'))
#   X <- X[is.element(it.X,it.y),]
#   nX <- sum(apply(cbind(y,X),1,FUN=function(x) x[1] < quantile(x[-1],probs=0.05) | x[1] > quantile(x[-1],probs=0.95)))
#   score <- pbinom(nX,size=length(y),prob=0.1)
#   return(c(score,lon(x),lat(x)))
# }

# 
# shinyServer(function(input, output, session) {
#   #countview <- reactiveValues(i = 1)
#   
#   ## Try to get the location names to depend on whether temperature of precipitation stations
# 
#     
#   ## Show map of gridded temperature
#   output$maps <- renderPlot({ 
#     it <- range(as.numeric(input$dates1))
#     is <- list(lon=as.numeric(input$lon1),lat=as.numeric(input$lat1))
#     season <- switch(tolower(as.character(input$season1)),
#                      'winter'=1,'spring'=2,'summer'=3,'autumn'=4)
#     rcp <- switch(tolower(as.character(input$rcp1)),
#                   'rcp4.5'=1,'rcp2.6'=2,'rcp8.5'=3)
#     param <- switch(tolower(as.character(input$param1)),
#                     'temperature'=0,'wet-day freq.'=12,'precip. intensity'=24,
#                     'precip. sum'=-1)
#     FUN <- switch(tolower(as.character(input$aspect)),
#                   "mean value"="mean", "change"="change","trend"="trend","variability"="sd")
#     FUNX <- switch(tolower(as.character(input$stats)),
#                    "ensemble mean"="mean", "ensemble spread"="sd")
#     
#     if (param>=0) li <- (rcp-1)*4+season + param else
#                   li <- (rcp-1)*4+season + 12
#     im <- is.element(gcmnames,input$im)
#     gcnames <<- names(Z4[[li]])[-c(1,2,length(Z4[[1]]))]
#     ## Check if thetotal precipitation needs to be estimated:
#     if (param>=0) {
#       zmap <- Z4[[li]]
#     } else {
#       ## The total precipitation = 90 * fw * mu which requires a transform from EOF to
#       ## a field object
#       zmap1 <- Z4[[(rcp-1)*4+season + 12]]
#       zmap2 <- Z4[[(rcp-1)*4+season + 24]]
#       itx <- it
#       if ((FUN=="change")) itx <- c(min(1961,itx),max(1990,itx)) else 
#       if (FUN=="trend") itx[2] <- max(c(itx[2],itx[1]+30))
#       zmap1 <- map(zmap1,FUN="mean",FUNX=FUNX,it=itx,is=is,im=im,plot=FALSE)
#       zmap2 <- map(zmap2,FUN="mean",FUNX=FUNX,it=itx,is=is,im=im,plot=FALSE)
#       zmap <- 90*zmap1*zmap2
#       zmap <- attrcp(zmap1,zmap)
#       class(zmap) <- class(zmap1)
#       attr(zmap,'variable') <- 'precip'
#       attr(zmap,'unit') <- 'mm/season'
#       rm('zmap1','zmap2')
#     }
#     
#     if (FUN=="trend") it[2] <- max(c(it[2],it[1]+30))
#     main <- paste('Downscaled',FUN,input$season1,tolower(input$param1),'for',it[1],'-',it[2],
#                   'following',toupper(input$rcp1),'based on',sum(im),'model runs')
#     if (FUN=="change") {
#       y <- map(zmap,FUN="mean",FUNX=FUNX,it=it,is=is,im=im,plot=FALSE)
#       it0 <- range(as.numeric(input$baseline))
#       #it0 <- c(1961,1990)
#       y0 <- map(zmap,FUN="mean",FUNX=FUNX,it=it0,is=is,im=im,plot=FALSE)
#       if (is.zoo(y)) coredata(y) <- t(coredata(t(y)) - apply(coredata(t(y0)),1,FUN='mean')) else
#                      y <- y - y0  
#     } else y <- map(zmap,FUN=FUN,FUNX=FUNX,it=it,is=is,im=im,plot=FALSE)
#     
#     map(y,main=main,FUN="mean",new=FALSE)
#     }, height=function(){0.8*session$clientData$output_maps_width})
#    
#  
#   ## Unfinished!  
#   output$map.quality <- renderPlot({ 
#     season <- switch(tolower(as.character(input$season6)),
#                      'winter'=1,'spring'=2,'summer'=3,'autumn'=4)
#     rcp <- switch(tolower(as.character(input$rcp6)),
#                   'rcp4.5'=1,'rcp2.6'=2,'rcp8.5'=3)
#     param <- switch(tolower(as.character(input$param6)),
#                     'temperature'=0,'wet-day freq.'=12,'precip. intensity'=24,
#                     'precip. sum'=-1)
#     li <- (rcp-1)*4+season+param
#     is <- list(lon=as.numeric(input$lon6),lat=as.numeric(input$lat6))
#     im <- is.element(gcmnames,input$im)
#     
#     ## If all models are selected, use pre-calculated quality scores else calculate for 
#     ## selected models
#     if (sum(im) < length(im)) {
#       zz <- Z4[[li]]
#       class(zz) <- c('dsensemble','pca','list')
#       ## Reduce the matrix size and pick one station before the recovery of the original format
#       im <- is.element(gcmnames,input$im)
#       zz <- subset(zz,is=is,im=im)
#       z <- as.station(zz)
#       if (as.character(input$quality6)=="trend") {
#         score <- unlist(lapply(z,trendscore)) 
#       } else if (as.character(input$quality6)=="spread") {
#         score<- unlist(lapply(z,varscore))
#       }
#       dim(score) <- c(3,length(z))
#       lons <- lon(zz$pca); lats <- lat(zz$pca)
#     } else {
#       ## Pre-calculated scores for whole ensemble
#       if (as.character(input$quality6)=="trend") {
#         score <- quality$trend[[li]]
#       } else {
#         score <- quality$range[[li]]
#       }
#       if (param==0) {
#         lons <- lon(quality$t2m); lats <- lat(quality$t2m)
#       } else {
#         lons <- lon(quality$pre); lats <- lat(quality$pre)
#       }
#       ixy <- (lons >= is$lon[1]) & (lons <= is$lon[2]) &
#              (lats >= is$lat[1]) & (lats <= is$lat[2])
#       score <- score[,ixy]; lons <- lons[ixy]; lats <- lats[ixy]
#     }
#     main <- paste('Quality:',input$quality6,' for ',
#                   input$season6,input$param6,' (',sum(im),input$rcp6,' runs)',sep='')
#     plot(lons,lats,xlab='',ylab='',main=main,lwd=2,cex=2,pch=19,col='grey80')
#     grid()
# 
#     data(geoborders)
#     lines(geoborders,col='grey')
#    
#     col <- rep(rgb(0,0.5,0.5),length(lons)); pch <- rep(19,length(lons))
#     q2 <- score[1,] < 0.1 | score[1,] > 0.9
#     q3 <- score[1,] < 0.05 | score[1,] > 0.95
#     col[q2] <- rgb(0.5,0.5,0); pch[q2] <- 15
#     col[q3] <- rgb(1,0,0);   pch[q3] <- 17
#     points(score[2,],score[3,],pch=pch,cex=1.5,col=col)
#     par0 <- par()
#     par(new=TRUE,fig=c(0.05,0.25,0.9,1),mar=rep(1,4),xaxt='n',yaxt='n')
#     image(cbind(1:3,1:3),col=c(rgb(0,0.5,0.5),rgb(0.5,0.5,0),rgb(1,0,0)))
#     text(0,1,'[10,90]',col='white',cex=0.8,pos=1)
#     text(0.50,1,'[05,95]',col='white',cex=0.8,pos=1)
#     text(1,1,'outside',col='white',cex=0.8,pos=1)
#     print("almost done")
#     par(par0)
#     },height=function(){0.8*session$clientData$output_map.quality_width} )#600})
#   
#   ## Show thedifference between one selected model and the mean of the rest of the ensemble.
#   output$plot1model <- renderPlot({ 
#     it <- range(as.numeric(input$dates7))
#     is <- list(lon=as.numeric(input$lon7),lat=as.numeric(input$lat7))
#     season <- switch(tolower(as.character(input$season7)),
#                      'winter'=1,'spring'=2,'summer'=3,'autumn'=4)
#     rcp <- switch(tolower(as.character(input$rcp7)),
#                   'rcp4.5'=1,'rcp2.6'=2,'rcp8.5'=3)
#     param <- switch(tolower(as.character(input$param7)),
#                     'temperature'=0,'wet-day freq.'=12,'precip. intensity'=24)
#     li <- (rcp-1)*4+season+param
#     gcnames <- names(Z4[[li]])[-c(1,2,length(Z4[[1]]))]
#     im1 <- is.element(gcmnames,input$im7)
#     im <-  !is.element(gcmnames,input$im7) & is.element(gcmnames,input$im) 
#     
#     z1 <- subset(Z4[[li]],im=im1,it=it,is=is)
#     zz <- subset(Z4[[li]],im=im,it=it,is=is)
#     
#     y1 <- map(z1,plot=FALSE)
#     yy <- map(zz,plot=FALSE)
#     coredata(y1) <- coredata(y1) - coredata(yy)
#     main <- paste(gcmnames[im1],' - ensemble mean (number of runs=',sum(im),') ',
#                   season,'/',input$rcp1,': ',it[1],'-',it[2],sep='')
#     map(y1,main=main,new=FALSE)
#     #plot(rnorm(100),main=main)
#    },height=function(){0.8*session$clientData$output_plot1model_width} )#600})
#   
# 
#   output$use.stats <- renderText({
#     txt <- paste(countview$i,"actions")
#   })
#   
# })

# names.gcms <- function(data="stats",project_id="CMIP5",var="tas") {
#   if(data=="stats") {
#     if(is.null(var)) var <- "tas"
#     nm <- names(stats[[var]]$ff)
#   } else if(data=="meta") {
#     nm <- apply(M,1,function(x) {
#       y <- paste(x["gcm"],x["rip"],sep=".")
#       if(!is.na(x["rcm"])) y <- paste(x["rcm"],y,sep=" / ")
#       return(y) } )
#     im <- rep(TRUE,length(nm))
#     if(!is.null(var)) im <- M$var==var
#     if(!is.null(project_id)) im <- M$project_id==project_id
#     im <- seq(length(nm))[im]
#     nm <- nm[im]
#     names(nm) <- im
#   } else if(data=="ceof") {
#     if(var %in% c("temp","tas")) {
#       nm <- ceof["tas"][project_id]
#     } else {
#       nm <- ceof["pr"][project_id]
#     }
#   }
#   return(nm)
# }

