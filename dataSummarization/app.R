# app.R
source('global.R')

ui <- fluidPage(
  dataSummarizationModuleUI("demo",projectDf,keyList)
)

server <- function(input,output,session){
  callModule(dataSummarizationModule,"demo",projectDf)
}

shinyApp(ui,server)