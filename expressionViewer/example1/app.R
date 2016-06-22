source("global.R")

ui <- fluidPage(
  expressionViewerModuleUI("demo") 
)

server <- function(input,output){
  callModule(expressionViewerModule,"demo",eset.data,"demo")
}

shinyApp(ui,server)