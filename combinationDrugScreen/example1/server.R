shinyServer(function(input, output, session) {
  callModule(combinationDrugScreenModule,id = "demo",session = session,combinedData,sampleInfo, tag = "demo")
})