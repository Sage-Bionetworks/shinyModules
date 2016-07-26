source("global.R")

ui <- fluidPage(
  expressionViewerModuleUI("demo") 
)

server <- function(input,output){
  callModule(expressionViewerModule,"demo",eset.data,pathways_list,"demo")
}

shinyApp(ui,server)