source("global.R")

ui <- fluidPage(
  heatmapModuleUI("demo") 
)

server <- function(input,output){
  callModule(heatmapModule,"demo",eset.data,"demo")
}

shinyApp(ui,server)