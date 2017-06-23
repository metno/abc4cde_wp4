
library(shiny)
if (!require("DT")) install.packages('DT')
library(DT)
library(DECM)

#data("statistics.cmip.tas.2071-2100")
#gcmnames <- paste(seq(length(store)),names(store),sep=". ")

navbarPage(title = 'Climate change scatterplot of temperature and precipitation', 
           id = 'x0',
           header=img(src="banner_c3s.png", 
                      style="width:700px; float:left; margin-bottom:50px; margin-left:20px"),
           #tabPanel('Scatterplot of climate change', 
                    #column(3,
                    #       checkboxGroupInput("gcms",
                    #                          label = "Climate models",
                    #                          choices = gcmnames)),
                    column(3,
                           selectInput("season",
                                       label = "Season",
                                       choices = c("Annual mean","Winter","Spring","Summer","Autumn"),
                                       selected = "Annual mean")),
                    column(3,
                           selectInput("period",
                                       label = "Future time period",
                                       choices = c("Far future (2071-2100)",
                                                   "Near future (2021-2050)"),
                                       selected = "Far future (2071-2100)")),
                    plotOutput("dtdpr")
           #)
  
)
