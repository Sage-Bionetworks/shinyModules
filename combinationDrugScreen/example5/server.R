shinyServer(function(input, output, session) {
  callModule(combinationDrugScreenModule,id = "demo",session = session,combinedData,tag = "demo")
})