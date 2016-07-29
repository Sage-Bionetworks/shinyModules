combinationDrugScreenModuleUI <- function(id,sampleInfo){
  ns <- NS(id)
  tagList(
    myHeader <- dashboardHeader(title="Combination Drug Screen", disable=TRUE),
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
            column(width = 3,
                   h4('1. Select Samples'),
                   selectInput(ns('selected_sample'),NULL, choices = unique(sampleInfo$sample),
                    selectize=T, multiple=F, selected = unique(sampleInfo$sample)[1])
            ),
             column(width = 2,
                   h4('2. Select format'),
                   selectInput(ns('selected_format'),NULL, choices = unique(sampleInfo$format),
                    selectize=T, multiple=F, selected = unique(sampleInfo$format)[1])
            ),
            column(width = 2,
                   h4('3. Select Metric'),
                   selectInput(ns('selected_metric'),NULL, choices = unique(sampleInfo$metric),
                    selectize=T, multiple=F, selected = unique(sampleInfo$metric)[1])
            ),
            column(width = 5,
                   h4('4. Select Drugs'),
                   selectInput(ns('selected_drug1'),"Drug 1", choices = NULL,selectize=T, multiple=F, selected = NULL),
                   selectInput(ns('selected_drug2'),"Drug 2", choices = NULL,selectize=T, multiple=F, selected = NULL)
            ),
            actionButton(ns("updateButton"), "Update")
        )
      ),
      
      fluidRow(
        box(width = 12,
            plotOutput(ns("plots"),height = "520px")
        )
      )
      )
    
      )
  dashboardPage(header=myHeader, sidebar=mySidebar, body=myBody,
                skin = "blue")
  } 

combinationDrugScreenModule <- function(input,output,session,combinedData,sampleInfo,tag){
  observeEvent(input$selected_sample,{
    sample <- input$selected_sample
    format <- unique(sampleInfo[sampleInfo$sample == sample,]$format)
    updateSelectInput(session,"selected_format",choices = format, selected = format[1])
  })

  observeEvent(input$selected_format,{
    sample <- input$selected_sample
    format <- input$selected_format
    metric <- unique(sampleInfo[sampleInfo$sample == sample & sampleInfo$format == format,]$metric)
    updateSelectInput(session,"selected_metric",choices = metric, selected = metric[1])

  })


  dataset <- reactive({
    sample <- input$selected_sample
    format <- input$selected_format
    metric <- input$selected_metric
    
    sampleName <- paste(sample,format,metric,sep = "_")
    dataset <- combinedData[[sampleName]]
    
    dataset
  })

sampleInfoDf <- reactive({
  dataset <- dataset()
  pairs <- dataset$drug.pairs
  #remove "null"
  pairs <- pairs[pairs$drug.row != "null" & pairs$drug.col != "null",]
  pairs
})

observeEvent(sampleInfoDf(),{
  pairs <- sampleInfoDf()
  drug.row <- sort(unique(pairs$drug.row))
  updateSelectInput(session,"selected_drug1",choices = drug.row, selected = drug.row[1])
})

observeEvent(input$selected_drug1,{
  pairs <- sampleInfoDf()
  drug.row <- input$selected_drug1
  drug.col <- sort(pairs[pairs$drug.row == drug.row,]$drug.col)
  updateSelectInput(session,"selected_drug2",choices = drug.col, selected = drug.col[1])
})
  
flt_dataset <- eventReactive(input$updateButton,{
    dataset <- dataset()
    drug1 <- input$selected_drug1
    drug2 <- input$selected_drug2
    
    pairs <- dataset$drug.pairs
    sampleDf <- pairs[pairs$drug.row == drug1 & pairs$drug.col == drug2,]
    blockMat <- dataset$dose.response.mats[[sampleDf$blockIDs]]
    
    list(sampleInfo = sampleDf, mat = blockMat)
  })
  
  
  output$plots <- renderPlot({
    validate(need(input$updateButton, "Please click \"Update\"."))
    dataset <- flt_dataset()
    sampleInfo <- dataset[["sampleInfo"]]
    response.mat <- dataset[["mat"]]
    
    num.row <- length(response.mat)
    data.plot <- data.frame(x = numeric(num.row), y = numeric(num.row),Inhibition = numeric(num.row))
    data.plot$Inhibition <- round(c(response.mat), 2)
    data.plot$y <- rep(c(1:ncol(response.mat)), nrow(response.mat))
    data.plot$x <- rep(1:nrow(response.mat), each = ncol(response.mat))
    data.plot$x <- as.factor(data.plot$x)
    data.plot$y <- as.factor(data.plot$y)
    conc.unit <- sampleInfo$concUnit ## concentration unit
    
    unit.text <- paste("(", conc.unit, ")", sep = "")
    
    drug.row <- sampleInfo$drug.row
    drug.col <- sampleInfo$drug.col
    
    plot.title <- "Heatmap"
    axis.x.text <- round(as.numeric(colnames(response.mat)), 1)
    axis.y.text <- round(as.numeric(rownames(response.mat)), 1)
    dose.response.p <- ggplot(data.plot, aes_string(x = "x", y = "y")) + geom_tile(aes_string(fill = 'Inhibition')) +
      geom_text(aes_string(fill = 'Inhibition', label = 'Inhibition')) +
      scale_fill_gradient2(low = "green", high = "red", midpoint = 0, name = "Inhibiton (%)") +
      scale_x_discrete(labels = axis.x.text) + scale_y_discrete(labels = axis.y.text) +
      xlab(paste(drug.col, unit.text, sep = " ")) + ylab(paste(drug.row, unit.text, sep = " "))
    dose.response.p <- dose.response.p + theme(axis.text.x = element_text(color = "red", face = "bold", size = 15))
    dose.response.p <- dose.response.p + theme(axis.text.y = element_text(color = "red", face = "bold", size = 15))
    dose.response.p <- dose.response.p + theme(axis.title = element_text(size=15))
    dose.response.p <- dose.response.p + ggtitle(plot.title)

    single.fitted <- FittingSingleDrug(response.mat)

    layout(matrix(c(1, 3, 2, 3), 2, 2, byrow = TRUE))
    # plot the curve for the row drug
    suppressWarnings(par(mgp=c(3, .5, 0)))
    x.lab <- paste("Concentration", unit.text, sep = " ")

    plot(single.fitted$drug.row.model, xlab = x.lab, ylab = "Inhibition (%)", type = "obs", col = "red", cex = 1.5, pch = 16, xtsty = "base5")
    plot(single.fitted$drug.row.model, xlab = x.lab, ylab = "Inhibition (%)", type = "none", cex = 1.5, add = T, lwd = 3)
    title(paste("Dose-response curve for:", drug.row), cex.main = 1)

    # plot the curve for the col drug
    plot(single.fitted$drug.col.model, xlab = x.lab, ylab = "Inhibition (%)", type = "obs", col = "red", cex = 1.5, pch = 16, xtsty = "base5")
    plot(single.fitted$drug.col.model, xlab = x.lab, ylab = "Inhibition (%)", type = "none", cex = 1.5, add = T, lwd = 3)
    title(paste("Dose-response curve for:", drug.col), cex.main = 1)

    plot.new()
    print(dose.response.p, vp = viewport(height = unit(1, "npc"), width = unit(0.5, "npc"), just = c("left","top"), y = 1, x = 0.5))
  })
 
}



