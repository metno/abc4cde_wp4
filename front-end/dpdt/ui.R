
library(shiny)
if (!require("DT")) install.packages('DT')
library(DT)
library(DECM)

#data("statistics.cmip.tas.2071-2100")
data("metaextract")
im <- meta$project_id=="CMIP5" & meta$var=="tas"
gcmnames <- paste(seq(sum(im)),": ",meta$gcm[im],".",meta$gcm_rip[im],sep="")

navbarPage(title = 'Projections of temperature and precipitation', 
           id = 'x0',
           header=img(src="banner_c3s.png", 
                      style="width:700px; float:left; margin-bottom:50px; margin-left:20px"),
           tabPanel('Scatterplot of climate change', 
                    column(3,
                           selectInput("season",
                                       label = "Season",
                                       choices = c("Annual mean","Winter","Spring","Summer","Autumn"),
                                       selected = "Annual mean"),
                           selectInput("period",
                                       label = "Future time period",
                                       choices = c("Far future (2071-2100)",
                                                   "Near future (2021-2050)"),
                                       selected = "Far future (2071-2100)"),
                           checkboxGroupInput("gcms",
                                              label = "Climate models (GCMs)",
                                              choices = gcmnames,
                                              selected = gcmnames)),
                    plotOutput("dtdpr")
           )
)
