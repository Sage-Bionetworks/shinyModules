source("global.R")

ui <- fluidPage(
  drugScreenModuleUI("demo") 
)

server <- function(input,output,session){
  #test <- callModule(testModule,"demo")
  callModule(drugScreenModule,id = "demo",session = session, summarizedData = summarizedData, tag = "demo")#,test)
}

shinyApp(ui,server)