shinyServer(function(input, output, session) {
  callModule(combinationDrugScreenModule,id = "demo",session = session,normData,sampleInfo, tag = "demo")
})