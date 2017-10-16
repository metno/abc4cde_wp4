## ABC4CDE/DEMC - R-script for prototype tool WP4
## andreas.dobler@met.no  Oslo, Norway, 2017-05-16
## 

##  RShiny server file

library(shiny)
library(plotrix)

#Load metadata from file
#AD: THIS IS CURRENTLY A LOCAL COPY. The final idea is to have it updated by the back-end
load("metaextracted.rda")

#Load standard deviations and correlations from file
#AD: THIS IS CURRENTLY A LOCAL COPY. The final idea is to have it updated by the back-end
load("sd_cor.rda")

shinyServer(function(input, output) {
		      #Source the taylor diagram function with correlation and standard devaition as input
		      source("taylor_sdR_in.R")
		        
		        #Taylor diagram
		        #AD: SHOWS CURRENTLY ALL MODELS WITH COR. AND SD INFORMATION!
		        output$taylor <- renderPlot({
				    
				    #Find models which have correlation and sd information:
				    #Get data URLs from metadata
				    # data <- data.frame(DataURL= sapply(mdList,function(X) gsub("dodsC","fileServer",X$dataurl)))
				    data <- data.frame(Project=sapply(mdList,function(X) X$globat$project_id),
						                              RCM=sapply(mdList,function(X) X$globat$model_id),
									                             GCM= sapply(mdList,function(X) X$globat$driving_model_id),
									                             DataLink= sapply(mdList,function(X) gsub("dodsC","fileServer",X$dataurl)),
												                            OPENDAP = sapply(mdList,function(X) paste0("<a href='",X$dataurl,".html ' target='_blank'>","OPENDAP access","</a>")))
				        
				        #Get data frame indices and put correlation, SSD and SSDRef into data frame
				        data_indx <- unlist(sapply(RSD_data$DataURL, function(X) which(gsub("dodsC","fileServer",X)  == data$DataLink)))
				        data$CorrWithRef[data_indx] <- RSD_data$CorrWithRef
					    data$SSD[data_indx] <- RSD_data$SSD
					    data$SSDRef[data_indx] <- RSD_data$SSDRef
					        
					        #Sub selection
					        if (input$project != "All") {
							      data <- data[data$Project == input$project,]
					        }
					        if (input$rcm != "All") {
							      data <- data[data$RCM == input$rcm,]
						    }
						    if (input$gcm != "All") {
							          data <- data[data$GCM == input$gcm,]
						        }
						    
						    #Create taylor plots
						    #Dummy plot
						    taylor.diagram_Rsd(1,1,2,col="white",main="Spatial field metrics",normalize = TRUE,ylab="Normalized standard devaition",pch=NA)
						        
						        #Add models
						        for (i in which(!is.na(data$CorrWithRef)))
								    {
									          taylor.diagram_Rsd(data$CorrWithRef[i],data$SSDRef[i],data$SSD[i],add=T,col=colors()[3*as.numeric(row.names(data))[i]],normalize=T)
						        }
						        
						        #add a legend
						        legend(2.8,2.8,legend=as.numeric(row.names(data))[which(!is.na(data$CorrWithRef))],pch=19,cex=2,col=colors()[3*as.numeric(row.names(data))[which(!is.na(data$CorrWithRef))]])
							    
							  })
		        
		        #Model selection table
		        output$table <- DT::renderDataTable(DT::datatable({
				    
				    #Get data selections and info from metadata
				    data <- data.frame(Project=sapply(mdList,function(X) X$globat$project_id),
						                              RCM=sapply(mdList,function(X) X$globat$model_id),
									                             GCM= sapply(mdList,function(X) X$globat$driving_model_id),
									                             DataLink= sapply(mdList,function(X) gsub("dodsC","fileServer",X$dataurl)),
												                            OPENDAP = sapply(mdList,function(X) paste0("<a href='",X$dataurl,".html ' target='_blank'>","OPENDAP access","</a>")))
				        
				        # data$DataURL = sapply(data$DataURL,function(X) paste0("<a href='",X,"' target='_blank'>",X,"</a>"))
				        
				        #Get data frame indices and put correlation, SSD and SSDRef into data frame
				        data_indx <- unlist(sapply(RSD_data$DataURL, function(X) which(gsub("dodsC","fileServer",X)  == data$DataLink)))
				        data$CorrWithRef[data_indx] <- RSD_data$CorrWithRef
					    data$SSD[data_indx] <- RSD_data$SSD
					    data$SSDRef[data_indx] <- RSD_data$SSDRef
					        
					        #Replace URL with hyperlink
					        data[,"DataLink"] <- sapply(data[,"DataLink"],function(X) paste0("<a href='",X,"' target='_blank'>","get data here","</a>"))
					        
					        #Sub selection
					        if (input$project != "All") {
							      data <- data[data$Project == input$project,]
						    }
						    if (input$rcm != "All") {
							          data <- data[data$RCM == input$rcm,]
						        }
						    if (input$gcm != "All") {
							          data <- data[data$GCM == input$gcm,]
						        }


						        #Display data
						        # data[,c("Project","RCM","GCM","CorrWithRef","SSD","SSDRef","DataLink")]#Incl. SSD, etc.
						        data[,c("Project","RCM","GCM","CorrWithRef","DataLink","OPENDAP")]
						      },options = list(lengthMenu = c(10, 20, 50, 100), pageLength = 20), escape = 1:4))
			  
			  
})


