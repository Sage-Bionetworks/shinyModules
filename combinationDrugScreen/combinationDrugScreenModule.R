combinationDrugScreenModuleUI <- function(id,combined_data){
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
            column(width = 4,
                   h4('1. Select Samples'),
                   selectInput(ns('selected_sample'),NULL, choices = unique(combined_data$sample),
                    selectize=T, multiple=T, selected = unique(combined_data$sample)[1:3])
            ),
            column(width = 3,
                   h4('2. Select Assay'),
                   selectInput(ns('selected_assay'),NULL, choices = NULL,
                    selectize=T, multiple=F, selected = NULL)
            ),
            column(width = 5,
                   h4('3. Select Drugs'),
                   selectInput(ns('selected_drug1'),"Drug 1", choices = NULL,selectize=T, multiple=F, selected = NULL),
                   selectInput(ns('selected_drug2'),"Drug 2", choices = NULL,selectize=T, multiple=F, selected = NULL)
            ),
            actionButton(ns("updateButton"), "Update")
        )
      ),
      
      fluidRow(
        tabBox(width = 12,height = "600px",
          tabPanel("Dose Response",
              uiOutput(ns("facet_by_ui")),
              plotOutput(ns("doseResp_plots"),height = "500px")
          ),
          tabPanel("Heatmap",
                plotOutput(ns("heatmap_plots"),height = "500px")
          )
        )
      )
      )
    
      )
  dashboardPage(header=myHeader, sidebar=mySidebar, body=myBody,
                skin = "blue")
} 

combinationDrugScreenModule <- function(input,output,session,combined_data,tag){
  ns <- NS(tag)
  observeEvent(input$selected_sample,{
    assay <- unique(filter(combined_data, sample %in% input$selected_sample)$assay)
    updateSelectInput(session,"selected_assay",choices = assay, selected = assay[1])
  })

  dataset <- reactive({
    dataset <- filter(combined_data, sample %in% input$selected_sample)
    dataset <- dataset[dataset$assay == input$selected_assay,]
    dataset
  })

  observeEvent(dataset(),{
    dataset <- dataset()
    drug.row <- sort(unique(dataset$drug1))
    updateSelectInput(session,"selected_drug1",choices = drug.row, selected = drug.row[1])
  })

  observeEvent(input$selected_drug1,{
    dataset <- dataset()
    drug.col <- sort(unique(dataset[dataset$drug1 == input$selected_drug1,]$drug2))
    updateSelectInput(session,"selected_drug2",choices = drug.col, selected = drug.col[1])
  })
    
  flt_dataset <- eventReactive(input$updateButton,{
    validate(need(!is.null(input$selected_sample), "At least one sample needs to be selected." ),
             need(length(input$selected_sample) <= 5, "You can select up to 5 samples." ))
    dataset <- dataset()

    flt_dataset <- dataset[dataset$drug1 == input$selected_drug1 & dataset$drug2 == input$selected_drug2,]
    flt_dataset
  })
  
  # Dose Response Tab
  output$facet_by_ui <- renderUI({
    if(input$updateButton){
      tagList(
        selectInput(ns("facet_by"),"Facet By", choices = c("sample","numDosagePoints"), selected = "numDosagePoints")
      )
    }
  })
  
  plot_dataset <- reactive({
    flt_dataset <- flt_dataset()
    
    # Get data for drug1
    drug1_dataset <- flt_dataset[flt_dataset$conc2 == 0,]
    drug1_dataset <- drug1_dataset[,-which(names(drug1_dataset) %in% c("drug2","conc2"))]
    names(drug1_dataset) <- sub("[12]","",names(drug1_dataset))

    # Get data for drug2
    drug2_dataset <- flt_dataset[flt_dataset$conc1 == 0,]
    drug2_dataset <- drug2_dataset[,-which(names(drug2_dataset) %in% c("drug1","conc1"))]
    names(drug2_dataset) <- sub("[12]","",names(drug2_dataset))
  
    # Combine and remove 0s
    drug_dataset <- rbind(drug1_dataset,drug2_dataset)
    drug_dataset$response <- as.numeric(drug_dataset$response)/100
    drug_dataset <- drug_dataset[drug_dataset$conc != 0,]
    
    drug_dataset
  })

  doseResp_data <- reactive({
    dt <- plot_dataset()
    var <- c('drug', 'sample', 'numDosagePoints')
    doseRespData <- plyr::ddply(.data=dt, .variables = var,.fun = tmp_iterator, .parallel = T)
    doseRespData$grp <- doseRespData$sample
    doseRespData
  })
  
  output$doseResp_plots <- renderPlot({
    validate(need(input$updateButton, "Please click \"Update\"."))
    drug_dt <- plot_dataset()
    drug_dt$numDosagePoints <- as.character(drug_dt$numDosagePoints)
    resp_dt <- doseResp_data()
    resp_dt$fittedY <- resp_dt$fittedY*100
    resp_dt$numDosagePoints <- as.character(resp_dt$numDosagePoints)
    
    facet_by <- input$facet_by
    color_by <- "sample"
    if(facet_by == color_by){
      color_by <- "numDosagePoints" 
    }
    
    rangeVal <- quantile(log10(drug_dt$conc))
    labelVal <- c(rangeVal[2], rangeVal[3], rangeVal[4])
    p <- ggplot(drug_dt, aes(x = log10(conc), y = response*100)) 
    p <- p + geom_point(aes_string(color=color_by)) 
    p <- p + scale_color_brewer(type = "qual", palette = 2, direction = 1)
    p <- p + geom_line(data = resp_dt, aes_string(x = "fittedX", y = "fittedY", colour = color_by, group = color_by))
    p <- p + facet_grid(as.formula(paste(facet_by,"~","drug"))) + theme_bw(base_size = 15)
    p <- p + geom_hline(aes(yintercept=50), color='grey50', linetype='dashed')
    p <- p + scale_x_continuous(breaks = labelVal, labels = sapply(labelVal, function(x) format(signif(10^x,digits = 2),scientific = T)))
    p <- p + xlab('conc (uM)') + ylab('cell viability %') 
  
    p
  })
  
  # Heatmap Tab
  heatmap_dataset <- reactive({
    flt_dataset <- flt_dataset()
    data_info <- flt_dataset[,c("sample","numDosagePoints")]
    data_info <- data_info[!duplicated(data_info),]
    row.names(data_info) <- apply(data_info,1,function(x){
      paste(x[1], str_trim(x[2]),sep="_")
    })
    
    mats <- apply(data_info,1,function(x){
      sample <- x[1]
      size <- str_trim(x[2])
      dt <- flt_dataset[flt_dataset$sample == sample & flt_dataset$numDosagePoints == size,]
      dt <- dt[,c("conc1","conc2","response")]
      dt <- cast(dt, conc1 ~ conc2, value = "response")
      row.names(dt) <- dt$conc1
      dt$conc1 <- NULL
      return(dt)
    })
    
    mats
  })

  heatmap_list <- reactive({
    plotlist <- list()
    heatmap_dataset <- heatmap_dataset()
    
    plotlist <- lapply(names(heatmap_dataset), function(x){
      dt <- as.matrix(heatmap_dataset[[x]])
      
      num.row <- length(dt)
      data.plot <- data.frame(x = numeric(num.row), y = numeric(num.row),
                              response = numeric(num.row))
      data.plot$response <- round(c(dt), 2)
      data.plot$y <- rep(c(1:ncol(dt)), nrow(dt))
      data.plot$x <- rep(1:nrow(dt), each = ncol(dt))
      data.plot$x <- as.factor(data.plot$x)
      data.plot$y <- as.factor(data.plot$y)
      data.info <- str_split(x,"_")[[1]]
      plot.title <- paste("Sample:",data.info[1],
                          "\nNumber of Dosage Points:",data.info[2])
      
      
      axis.x.text <- signif(as.numeric(colnames(dt)), 2)
      axis.y.text <- signif(as.numeric(rownames(dt)), 2)
      p <- ggplot(data.plot, aes_string(x = "x", y = "y")) + geom_tile(aes_string(fill = 'response')) + 
        geom_text(aes_string(fill = 'response', label = 'response')) +
        scale_fill_gradient2(low = "green", high = "red", midpoint = 0, name = "response (%)") +
        scale_x_discrete(labels = axis.x.text) + scale_y_discrete(labels = axis.y.text) +
        xlab(paste(input$selected_drug1, "(uM)", sep = " ")) + ylab(paste(input$selected_drug2, "(uM)", sep = " "))
      
      p <- p + theme(axis.text.x = element_text(color = "red", face = "bold", size = 12))
      p <- p + theme(axis.text.y = element_text(color = "red", face = "bold", size = 12))
      p <- p + theme(axis.title = element_text(size=15))
      p <- p + ggtitle(plot.title)
      
      p
    })
    
    do.call(grid.arrange, c(plotlist, list(ncol = 2)))
  })
  

  output$heatmap_plots <- renderPlot({
    validate(need(input$updateButton, "Please click \"Update\"."))
    heatmap_list()
  })

}

get_drugResponse_stats <- function(conc,viability,...){
  res <- nplr(conc, viability,...)
  results <- getAUC(res)
  ICx_est = getEstimates(res, targets= c(.10,.20,.30,.40,.50,.60,.70,.80,.90))
  results['IC50'] = ICx_est[5,'x']
  results['maxResp'] = max(getYcurve(res)) #get the maximum efficacy of the drug
  
  fittedVals <- data.frame(fittedX = getXcurve(res),
                           fittedY = getYcurve(res))
  results <- cbind(results,fittedVals)
  results
}

tmp_iterator <- function(df){
  tryCatch({
    stats <- get_drugResponse_stats(df$conc, df$response, useLog=T)  
  },error=function(e){
    print(dim(df))
    print(df$conc)
    print(df$response)
    print(unique(df$sample))
    print(unique(df$drug))
    stop('stopped')
  })
}