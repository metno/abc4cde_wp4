
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
                      style="width:800px; float:left; margin-bottom:50px; margin-left:20px"),
           tabPanel('Seasonal cycles', 
                    column(3,
                           selectInput("variable",
                                       label = "Meteorological variable",
                                       choices = c("Temperature",
                                                   "Precipitation"),
                                       selected = "Temperature"),
                           selectInput("period",
                                       label = "Time period",
                                       choices = c("Present day (1981-2010)",
                                                   "Far future (2071-2100)",
                                                   "Near future (2021-2050)"),
                                       selected = "Present day (1981-2010)"),
                           checkboxGroupInput("gcms",
                                              label = "Climate models (GCMs)",
                                              choices = gcmnames,
                                              selected = gcmnames[[1]]),
                           numericInput("y1","y-axis max", NULL, min = -50, max = 50),#,
                           numericInput("y0","y-axis min", NULL, min = -50, max = 50)),
                           #sliderInput("ylim", "range of y-axis:",
                           #           min = -50, max = 50, value = c(NULL,NULL), step=0.2)),
                    column(6,
                           plotOutput("seasoncycle")),
                    column(3,
                           selectInput("region",
                                label = "Region",
                                choices = c("Global",
                                            "Alaska/N.W. Canada [ALA:1]","Amazon [AMZ:7]",
                                            "Central America/Mexico [CAM:6]","small islands regions Caribbean",
                                            "Central Asia [CAS:20]","Central Europe [CEU:12]",
                                            "Canada/Greenland/Iceland [CGI:2]","Central North America [CNA:4]",
                                            "East Africa [EAF:16]","East Asia [EAS:22]",
                                            "East North America [ENA:5]","South Europe/Mediterranean [MED:13]",
                                            "North Asia [NAS:18]","North Australia [NAU:25]",
                                            "North-East Brazil [NEB:8]","North Europe [NEU:11]",
                                            "Southern Africa [SAF:17]","Sahara [SAH:14]",
                                            "South Asia [SAS:23]","South Australia/New Zealand [SAU:26]",
                                            "Southeast Asia [SEA:24]","Southeastern South America [SSA:10]",
                                            "Tibetan Plateau [TIB:21]","West Africa [WAF:15]",
                                            "West Asia [WAS:19]","West North America [WNA:3]",
                                            "West Coast South America [WSA:9]","Antarctica",
                                            "Arctic","Pacific Islands region[2]",
                                            "Southern Topical Pacific","Pacific Islands region[3]",
                                            "West Indian Ocean")),
                                plotOutput("map")
			)
           )
)
