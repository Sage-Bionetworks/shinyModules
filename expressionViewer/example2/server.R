shinyServer(function(input,output){
  callModule(expressionViewerModule,"demo",eset.data,pathways_list,"demo")
})
