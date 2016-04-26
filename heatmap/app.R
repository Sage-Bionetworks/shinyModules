source("global.R")

ui <- fluidPage(
  heatmapModuleUI("test") 
)

server <- function(input,output){
  callModule(heatmapModule,"test",eset.data)
}

shinyApp(ui,server)