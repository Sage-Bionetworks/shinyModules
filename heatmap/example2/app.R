source("global.R")

ui <- fluidPage(
  heatmapModuleUI("demo") 
)

server <- function(input,output){
  callModule(heatmapModule,"demo",eset.data)
}

shinyApp(ui,server)