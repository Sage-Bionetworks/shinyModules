source("global.R")

ui <- fluidPage(
  drugScreenModuleUI("demo") 
)

server <- function(input,output){
  #test <- callModule(testModule,"demo")
  callModule(drugScreenModule,"demo",drugScreenData)#,test)
}

shinyApp(ui,server)