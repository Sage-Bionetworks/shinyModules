source("global.R")

ui <- fluidPage(
  drugScreenModuleUI(id = "demo", data = summarizedData) 
)

server <- function(input,output,session){
  #test <- callModule(testModule,"demo")
  callModule(drugScreenModule,id = "demo",session = session, summarizedData = summarizedData, tag = "demo")#,test)
}

shinyApp(ui,server)
