source("global.R")

ui <- fluidPage(
  drugScreenModuleUI("demo") 
)

server <- function(input,output){
  callModule(drugScreenModule,"demo",data1 = drug_normViab, data2 = drug_ICVals)
}

shinyApp(ui,server)