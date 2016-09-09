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

                                .col-sm-6:nth-child(odd){
                                  padding-bottom: 0;
                                  padding-top: 0;
                                }

                                .col-sm-6:nth-child(even){
                                  padding-top: 0;
                                }
                                '))),
      fluidRow(
        box(width=12, status='primary', collapsible=TRUE, 
            collapsed=FALSE, solidHeader=TRUE,
            title = tagList(shiny::icon("th-list", lib="glyphicon"),
                            "Data Selection"),
            column(width = 4,
                   h4('1. Select Samples'),
                   selectInput(ns('selected_sample'),NULL, choices = unique(sampleInfo$sample),
                    selectize=T, multiple=F, selected = unique(sampleInfo$sample)[1])
            ),
            column(width = 3,
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
        box(width = 6,
            plotOutput(ns("doseResp_plot1"),height = "520px")
            #helpText("Dose response plots")
        ),
        box(width = 6,
            plotOutput(ns("doseResp_plot2"),height = "520px")
        ),
        box(width = 6,
            plotOutput(ns("heatmap_plot1"),height = "520px")
        ),
        box(width = 6,
            plotOutput(ns("heatmap_plot2"),height = "520px")
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
    metric <- unique(sampleInfo[sampleInfo$sample == sample,]$metric)
    updateSelectInput(session,"selected_metric",choices = metric, selected = metric[1])
  })

  dataset <- reactive({
    sample <- input$selected_sample
    metric <- input$selected_metric
    
    format <- unique(sampleInfo$format)
    
    sampleName1 <- paste(sample,format[1],metric,sep = "_")
    sampleName2 <- paste(sample,format[2],metric,sep = "_")
    dataset <- combinedData[c(sampleName1,sampleName2)]
    dataset <- dataset[!is.na(names(dataset))]
    
    dataset
  })

  sampleInfoDf <- reactive({
    dataset <- dataset()

    pairs <- lapply(c(1:length(dataset)), function(x){
      tempName <- names(dataset)[x]
      temp <- dataset[[x]]$drug.pairs
      temp <- temp[temp$drug.row != "null" & temp$drug.col != "null",]
      temp$sampleName <- tempName
      return(temp)
    })

    if(length(pairs) == 2){
      pairs <- merge(pairs[[1]], pairs[[2]], all = TRUE)
    }else{
      pairs <- pairs[[1]]
    }

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
    drug.col <- sort(unique(pairs[pairs$drug.row == drug.row,]$drug.col))
    updateSelectInput(session,"selected_drug2",choices = drug.col, selected = drug.col[1])
  })
    
  flt_dataset <- eventReactive(input$updateButton,{
    dataset <- dataset()
    drug1 <- input$selected_drug1
    drug2 <- input$selected_drug2
    
    pairs <- sampleInfoDf()

    sampleDf <- pairs[pairs$drug.row == drug1 & pairs$drug.col == drug2,]

    result <- lapply(sampleDf$sampleName,function(x){
      temp <- dataset[[x]]
      tempInfoDf <- sampleDf[sampleDf$sampleName == x,]
      blockMat <- temp$dose.response.mats[[tempInfoDf$blockIDs]]
      return(list(sampleInfo = tempInfoDf,mat = blockMat))
    })

    result
  })
    
  output$doseResp_plot1 <- renderPlot({
    validate(need(input$updateButton, "Please click \"Update\"."))
    dataset1 <- flt_dataset()[[1]]
    sampleInfo1 <- dataset1[["sampleInfo"]]
    response.mat1 <- dataset1[["mat"]]
    single.fitted1 <- FittingSingleDrug(response.mat1)
    
    if(length(flt_dataset()) == 2){
      dataset2 <- flt_dataset()[[2]]
      sampleInfo2 <- dataset2[["sampleInfo"]]
      response.mat2 <- dataset2[["mat"]]
      single.fitted2 <- FittingSingleDrug(response.mat2)
    }
    
    conc.unit <- sampleInfo1$concUnit ## concentration unit
    unit.text <- paste("(", conc.unit, ")", sep = "")
    
    drug.row <- sampleInfo1$drug.row
    drug.col <- sampleInfo1$drug.col
    
    #layout(matrix(c(1, 3, 2, 3), 2, 2, byrow = TRUE))
    layout(matrix(c(1, 2), 2, 1, byrow = TRUE))
    x.lab <- paste("Concentration", unit.text, sep = " ")

    plot(single.fitted1$drug.row.model, xlab = x.lab, ylab = "Inhibition (%)", type = "obs", col = "red", cex = 1.5, pch = 16, xtsty = "base5")
    plot(single.fitted1$drug.row.model, xlab = x.lab, ylab = "Inhibition (%)", type = "none", cex = 1.5, add = T, lwd = 3)
    title(paste("Dose-response curve for:", drug.row), cex.main = 1)

    # plot the curve for the col drug
    plot(single.fitted1$drug.col.model, xlab = x.lab, ylab = "Inhibition (%)", type = "obs", col = "red", cex = 1.5, pch = 16, xtsty = "base5")
    plot(single.fitted1$drug.col.model, xlab = x.lab, ylab = "Inhibition (%)", type = "none", cex = 1.5, add = T, lwd = 3)
    title(paste("Dose-response curve for:", drug.col), cex.main = 1)
    
  })
  
  output$doseResp_plot2 <- renderPlot({
    validate(need(input$updateButton, "."))
    validate(need(length(flt_dataset()) == 2, "The second plot is not available for selected data."))
    
    dataset2 <- flt_dataset()[[2]]
    sampleInfo2 <- dataset2[["sampleInfo"]]
    response.mat2 <- dataset2[["mat"]]
    single.fitted2 <- FittingSingleDrug(response.mat2)
    
    conc.unit <- sampleInfo2$concUnit ## concentration unit
    unit.text <- paste("(", conc.unit, ")", sep = "")
    
    drug.row <- sampleInfo2$drug.row
    drug.col <- sampleInfo2$drug.col
    
    #layout(matrix(c(1, 3, 2, 3), 2, 2, byrow = TRUE))
    layout(matrix(c(1, 2), 2, 1, byrow = TRUE))
    x.lab <- paste("Concentration", unit.text, sep = " ")

    plot(single.fitted2$drug.row.model, xlab = x.lab, ylab = "Inhibition (%)", type = "obs", col = "red", cex = 1.5, pch = 16, xtsty = "base5")
    plot(single.fitted2$drug.row.model, xlab = x.lab, ylab = "Inhibition (%)", type = "none", cex = 1.5, add = T, lwd = 3)
    title(paste("Dose-response curve for:", drug.row), cex.main = 1)

    # plot the curve for the col drug
    plot(single.fitted2$drug.col.model, xlab = x.lab, ylab = "Inhibition (%)", type = "obs", col = "red", cex = 1.5, pch = 16, xtsty = "base5")
    plot(single.fitted2$drug.col.model, xlab = x.lab, ylab = "Inhibition (%)", type = "none", cex = 1.5, add = T, lwd = 3)
    title(paste("Dose-response curve for:", drug.col), cex.main = 1)
  })

  output$heatmap_plot1 <- renderPlot({
    validate(need(input$updateButton, "."))
    dataset <- flt_dataset()[[1]]
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
    dose.response.p
  })

  output$heatmap_plot2 <- renderPlot({
    validate(need(input$updateButton, "."))
    #validate(need(length(flt_dataset()) == 2, "The second heatmap is not available for selected data."))
    validate(need(length(flt_dataset()) == 2, "."))
    dataset <- flt_dataset()[[2]]
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
    dose.response.p
  })
  
#   output$plots <- renderPlot({
#     validate(need(input$updateButton, "Please click \"Update\"."))
#     dataset <- flt_dataset()
#     sampleInfo <- dataset[["sampleInfo"]]
#     response.mat <- dataset[["mat"]]
#     
#     num.row <- length(response.mat)
#     data.plot <- data.frame(x = numeric(num.row), y = numeric(num.row),Inhibition = numeric(num.row))
#     data.plot$Inhibition <- round(c(response.mat), 2)
#     data.plot$y <- rep(c(1:ncol(response.mat)), nrow(response.mat))
#     data.plot$x <- rep(1:nrow(response.mat), each = ncol(response.mat))
#     data.plot$x <- as.factor(data.plot$x)
#     data.plot$y <- as.factor(data.plot$y)
#     conc.unit <- sampleInfo$concUnit ## concentration unit
#     
#     unit.text <- paste("(", conc.unit, ")", sep = "")
#     
#     drug.row <- sampleInfo$drug.row
#     drug.col <- sampleInfo$drug.col
#     
#     plot.title <- "Heatmap"
#     axis.x.text <- round(as.numeric(colnames(response.mat)), 1)
#     axis.y.text <- round(as.numeric(rownames(response.mat)), 1)
#     dose.response.p <- ggplot(data.plot, aes_string(x = "x", y = "y")) + geom_tile(aes_string(fill = 'Inhibition')) +
#       geom_text(aes_string(fill = 'Inhibition', label = 'Inhibition')) +
#       scale_fill_gradient2(low = "green", high = "red", midpoint = 0, name = "Inhibiton (%)") +
#       scale_x_discrete(labels = axis.x.text) + scale_y_discrete(labels = axis.y.text) +
#       xlab(paste(drug.col, unit.text, sep = " ")) + ylab(paste(drug.row, unit.text, sep = " "))
#     dose.response.p <- dose.response.p + theme(axis.text.x = element_text(color = "red", face = "bold", size = 15))
#     dose.response.p <- dose.response.p + theme(axis.text.y = element_text(color = "red", face = "bold", size = 15))
#     dose.response.p <- dose.response.p + theme(axis.title = element_text(size=15))
#     dose.response.p <- dose.response.p + ggtitle(plot.title)
# 
#     single.fitted <- FittingSingleDrug(response.mat)
# 
#     layout(matrix(c(1, 3, 2, 3), 2, 2, byrow = TRUE))
#     # plot the curve for the row drug
#     suppressWarnings(par(mgp=c(3, .5, 0)))
#     x.lab <- paste("Concentration", unit.text, sep = " ")
# 
#     plot(single.fitted$drug.row.model, xlab = x.lab, ylab = "Inhibition (%)", type = "obs", col = "red", cex = 1.5, pch = 16, xtsty = "base5")
#     plot(single.fitted$drug.row.model, xlab = x.lab, ylab = "Inhibition (%)", type = "none", cex = 1.5, add = T, lwd = 3)
#     title(paste("Dose-response curve for:", drug.row), cex.main = 1)
# 
#     # plot the curve for the col drug
#     plot(single.fitted$drug.col.model, xlab = x.lab, ylab = "Inhibition (%)", type = "obs", col = "red", cex = 1.5, pch = 16, xtsty = "base5")
#     plot(single.fitted$drug.col.model, xlab = x.lab, ylab = "Inhibition (%)", type = "none", cex = 1.5, add = T, lwd = 3)
#     title(paste("Dose-response curve for:", drug.col), cex.main = 1)
# 
#     plot.new()
#     print(dose.response.p, vp = viewport(height = unit(1, "npc"), width = unit(0.5, "npc"), just = c("left","top"), y = 1, x = 0.5))
#   })

}