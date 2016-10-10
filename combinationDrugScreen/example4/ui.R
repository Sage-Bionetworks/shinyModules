#source("global.R")

shinyUI(fluidPage(
  combinationDrugScreenModuleUI(id = "demo",combinedData) 
))