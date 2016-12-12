dataSummarizationModuleUI <- function(id,projectDf,keyList){
  ns <- NS(id)
  tagList(
    myHeader <- dashboardHeader(title="Data Summarization", disable=TRUE),
    mySidebar <- dashboardSidebar(disable=TRUE),
    myBody <- dashboardBody(
      tags$head(tags$style(HTML('
                                .btn {
                                float:right;
                                }
                                '))),
      fluidRow(
        box(width=12, status='primary', collapsible=TRUE, 
            collapsed=FALSE, solidHeader=TRUE,
            title = tagList(shiny::icon("th-list", lib="glyphicon"),
                            "Data Selection"),
            column(width = 4,
                   h4('1. Project Scope'), 
                   selectInput(ns('selected_projects'),NULL, choices = projectDf$projectName,
                               selectize=T, multiple=T, selected = NULL)
            ),
            column(width = 4,
                   h4('2. Descriptive Keys'),
                   selectInput(ns('selected_des_keys'),NULL, choices = keyList,
                               selectize=T, multiple=T, selected = NULL)
            ),
            column(width = 4,
                   h4('3. Summary Keys'),
                   selectInput(ns('selected_sum_keys'),NULL, choices = keyList,
                               selectize=T, multiple=T, selected = NULL)
            ),
            actionButton(ns("generate_button"), "Generate")
        )
      ),
      
      fluidRow(
        tabBox(width = 12,
               tabPanel("Result",
                        downloadButton(ns("download_result")),
                        br(),
                        br(),
                        dataTableOutput(ns("summarized_table"))
               )
        )
      )
      )
  )

  dashboardPage(header=myHeader, sidebar=mySidebar, body=myBody,
                skin = "blue")
}


dataSummarizationModule <- function(input,output,session,projectDf){
  combinedTable <- reactive({
    annotationTbls <- projectDf[projectDf$projectName %in% input$selected_projects,]$annotationTable
    combinedTableList <- lapply(annotationTbls, function(x){
      tbl <- synTableQuery(paste("SELECT * FROM",x))
      tbl@values
    })
    combinedTable <- rbind.fill(combinedTableList)
    combinedTable
  })
  
  observe({
    keys <- colnames(combinedTable())
    updateSelectInput(session,'selected_des_keys',choices = keys)
  })
  
  observeEvent(input$selected_des_keys,{
    keys <- colnames(combinedTable())
    desKeys <- input$selected_des_keys
  
    updateSelectInput(session,"selected_sum_keys",choices = keys[! keys %in% desKeys])
  })
  
  summarizedTable <- eventReactive(input$generate_button,{
    combinedTable <- combinedTable()
    
    desKeys <- input$selected_des_keys
    sumKeys <- input$selected_sum_keys
    
    tbl <- ddply(combinedTable,desKeys,function(x){
      res <- sapply(sumKeys, function(key){
        length(unique(x[,key]))
      })
      res
    })
    
    tbl
  })
  
  output$summarized_table <- renderDataTable({
    summarizedTable()
  })
  
    output$download_result <- downloadHandler(
      filename = function() { paste(Sys.Date(),'summarizedTable.csv',sep="_") },
      content = function(file) {
        write.csv(summarizedTable(), file,row.names = FALSE)
      }
    )
}




