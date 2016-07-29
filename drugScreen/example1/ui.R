#source("global.R")

shinyUI(fluidPage(
  drugScreenModuleUI(id = "demo", data = summarizedData) 
 )
)