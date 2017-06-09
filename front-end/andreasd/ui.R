## ABC4CDE/DEMC - R-script for prototype tool WP4
## andreas.dobler@met.no  Oslo, Norway, 2017-05-16
## 

##  RShiny user intereface file
library(shiny)

#Load metadata from file
#AD: THIS IS CURRENTLY A LOCAL COPY. The final idea is to have it updated by the back-end
load("metaextracted.rda")

#Extract lists for selections from metadata
project_ids <- unique(sapply(mdList,function(X) X$globat$project_id))
model_ids <- unique(sapply(mdList,function(X) X$globat$model_id))
driving_model_ids <- unique(sapply(mdList,function(X) X$globat$driving_model_id))

#Make a simple data selector
##AD: NEED TO ADD VARIABLE SELECTION HERE!
navbarPage("Data selector",
  header=img(src="banner_c3s.png", 
   style="width:700px; float:left; margin-bottom:50px; margin-left:20px"),
  tabPanel("This is a prototype, e.g. no distinction btw. variables so far.", 
  #Data selection drop-downs
    fluidRow(column(3,selectInput("project",label = "Project",c("All",project_ids))),
             column(3,selectInput("rcm",label = "RCM",c("All",model_ids))),       
             column(3,selectInput("gcm",label = "GCM",c("All",driving_model_ids))),
      # Create a new row to display a taylor diagram and table
      fluidRow(column(4,plotOutput("taylor",height = "600px")),
               column(6,DT::dataTableOutput("table")))
    )
  )
)