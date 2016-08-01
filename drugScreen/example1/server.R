shinyServer(function(input,output,session){
  #test <- callModule(testModule,"demo")
  callModule(drugScreenModule,id = "demo",session = session, summarizedData = summarizedData, tag = "demo")#,test)
})
