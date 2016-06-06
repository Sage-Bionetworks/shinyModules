drugScreenModuleUI <- function(id){
  ns <- NS(id)
  tagList(
    myHeader <- dashboardHeader(title="Drug Screen Demo", disable=TRUE),
    mySidebar <- dashboardSidebar(disable=TRUE),
    myBody <- dashboardBody(
      fluidRow(
        box(width=12, status='primary', collapsible=TRUE, 
            collapsed=FALSE, solidHeader=TRUE,
            title = tagList(shiny::icon("th-list", lib="glyphicon"),
                            "Data Selection"),
            #tags$button(tags$style(".rightAlign{clear:both; float:right;}")),
            column(width = 5,
              #textOutput(ns("testTxt")),
              h4('1. Select Samples'), 
              selectInput(ns('samples'),NULL, choices = c('ALL', unique(drugScreenData$sample)),
                          selectize=T, multiple=T)
            ),
            column(width = 7,
              h4('2. Select Drugs'),
              h5("Select by drug name"),
              selectInput(ns('selected_drugs'),NULL, choices = c('ALL', unique(drugScreenData$drug)),
                          selectize=T, multiple=T),
              h5("Select by target class"),
              selectInput(ns('selected_class'),NULL, choices = c('ALL', unique(drugScreenData$target)),
                          selectize=T, multiple=T)
            ),
            actionButton(ns("updateButton"), "Update")
        ),
        box(width = 12, status = 'warning', collapsible = TRUE,
            collapsed = TRUE, solidHeader = TRUE,
            title = tagList(shiny::icon("filter", lib="glyphicon"),
                            "Filter"),
            column(width=3,
              #conditionalPanel(
              #  condition = "output.ic50Data == true",
                # IC50 filter
                sliderInput(ns('ic50_filter'), 'IC50',
                          min=10, max=100, value=c(1, 1000), 
                          step=10, round=0)
              #)
            ),
            column(width=3,
              # AC50 filter
              sliderInput(ns('ac50_filter'), 'AC50', 
                         min=10, max=100, value=c(1, 1000), 
                         step=10, round=0)
            ),
            column(width=3,
              # Curve class filter
              sliderInput(ns('curveClass'), 'Curve Class', 
                          min=-10, max=10, value=sort(unique(drugScreenData$curveClass)), 
                          step=1, round=0)
            ),
            column(width=3,
               # Max Response filter 
               sliderInput(ns('maxResp'), 'Max Response', 
                           min=1, max=100, value=c(1, 100), 
                           step=10, round=0)
            )
        )
      ),
             
      fluidRow(
             tabBox(width = 12,
               tabPanel("Max Response",
                 plotOutput(ns("drug_max_resp"),height="700px",width="auto",hoverId=NULL)#,
                 #helpText("ps: Efficacy is the % of cells that were killed by the given drug")
              ),
              tabPanel("IC50",
                 plotOutput(ns("drugScreen_IC50_plot"),height="700px",width="auto",hoverId=NULL)
              ),
              tabPanel("AC50",
                 plotOutput(ns("drugScreen_AC50_plot"),height="700px",width="auto",hoverId=NULL)
              ),
              tabPanel("QC",
                h5("Density Histograms Plots"),
                uiOutput(ns("QC_plots"))
#                 column(width = 6,
#                   h6("AUC"),
#                   plotOutput(ns("drugScreen_QC_plot1"),height="500px",width="auto",hoverId=NULL),
#                   h6("Max Response"),
#                   plotOutput(ns("drugScreen_QC_plot2"),height="500px",width="auto",hoverId=NULL)
#                 ),
#                 column(width = 6,
#                   h6("AC50"),
#                   plotOutput(ns("drugScreen_QC_plot3"),height="500px",width="auto",hoverId=NULL),
#                   conditionalPanel(
#                     condition = "output.icData",
#                     h6("IC50"),
#                     plotOutput(ns("drugScreen_QC_plot4"),height="500px",width="auto",hoverId=NULL)
#                   )
#                 )
              ),
              tabPanel("Data",
                downloadButton(ns("downloadData")),
                dataTableOutput(ns("drugScreen_dataTable"))
              )
             )
      )
    )
  )
  dashboardPage(header=myHeader, sidebar=mySidebar, body=myBody,
                skin = "blue")
} 

drugScreenModule <- function(input,output,session,data){
#   output$testTxt <- reactive({
#     test()
#   })
  
  get_selected_samples <- reactive({
    samples <- if('ALL' %in% input$samples) unique(data$sample) else input$samples
    validate(need(length(samples) != 0, "At least one cellLine needs to be selected" ) )
    samples
  })
  
  get_selected_drugs <- reactive({
    drugNames <- if('ALL' %in% input$selected_drugs) unique(data$drug) else input$selected_drugs
    targetClass <- if('ALL' %in% input$selected_class) unique(data$drug) else unique(data[data$target ==input$selected_class,]$drug)
    drugs <- union(drugNames,targetClass)
    validate(need(length(drugs) != 0, "At least one drug or target class needs to be selected"))
    drugs
  })
  
  get_drug_data <- eventReactive(input$updateButton,{
      flt_drug_data <- filter(data, drug %in% get_selected_drugs())  
      flt_drug_data <- filter(flt_drug_data, sample %in% get_selected_samples())  
      return(flt_drug_data)
  })
  
#   output$ic50Data <- reactive({
#     "IC50" %in% colnames(data)
#   })
#   
#   output$ac50Data <- reactive({
#     "AC50" %in% colnames(data)
#   })
  
  observe({
    drug_data <- get_drug_data()
    # AC50 slider
    ac50_min <- floor(min(drug_data$AC50,na.rm = T)*1000)/1000
    ac50_max <- floor(max(drug_data$AC50,na.rm = T)*1000)/1000
    updateSliderInput(session, "ac50_filter", max = ac50_max, value = c(10,ac50_max))
    
     # IC50 slider
     #icc50_min <- floor(min(drug_data$IC50,na.rm = T)*1000)/1000
     #ic50_max <- floor(max(drug_data$IC50,na.rm = T)*1000)/1000
     #updateSliderInput(session, "ic50_filter", min = ic50_min, max = ic50_max)
  })

  get_filtered_drug_data <- reactive({
    drug_data <- get_drug_data()
    
  })
  
  output$drug_max_resp <- renderPlot({
    validate(need(input$updateButton, "Please select samples and drugs, then click \"Update\" button." ))
    flt_drug_data <- get_drug_data()
    drug_levels <- flt_drug_data %>%
      group_by(drug) %>%
      summarise(med=median(maxResp, na.rm=T)) %>%
      arrange(desc(med)) %>% select(drug)
    drug_levels <- drug_levels$drug
    flt_drug_data$drug <- factor(flt_drug_data$drug,levels=drug_levels)
    p <- ggplot(data=flt_drug_data, aes(x=drug, y=maxResp, group=sample)) 
    p <- p + geom_point(aes(color=sample), size=3) + theme_bw()
    p <- p + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab('Drug') + ylab('Response')
    p
  })
  
  output$drugScreen_IC50_plot <- renderPlot({
    validate(need(!all(is.na(data$IC50)), "IC50 data does not exist." ))
    flt_drug_data <- get_drug_data()
    #remove NA
    flt_drug_data <- flt_drug_data[! is.na(flt_drug_data["IC50"]), ]
    #convert to log10
    flt_drug_data["IC50"] <- log10(as.numeric(flt_drug_data[,"IC50"])) 
    
    drug_levels <- flt_drug_data %>%
      group_by(drug) %>%
      summarise(med=median("IC50", na.rm=T)) %>%
      arrange(desc(med)) %>% select(drug)
    drug_levels <- drug_levels$drug
    flt_drug_data$drug <- factor(flt_drug_data$drug,levels=drug_levels)
    p <- ggplot(data=flt_drug_data, aes_string(x="drug", y=IC50, group="sample")) 
    p <- p + geom_point(aes(color=cellLine), size=3) + theme_bw()
    p + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab('Drug') + ylab('IC50')
  })
  
  output$drugScreen_AC50_plot <- renderPlot({
    validate(need(!all(is.na(data$AC50)), "AC50 data does not exist." ))
    flt_drug_data <- get_drug_data()
    #remove NA
    flt_drug_data <- flt_drug_data[! is.na(flt_drug_data["AC50"]), ]
    #convert to log10
    flt_drug_data["AC50"] <- log10(as.numeric(flt_drug_data[,"AC50"])) 
    
    drug_levels <- flt_drug_data %>%
      group_by(drug) %>%
      summarise(med=median("AC50", na.rm=T)) %>%
      arrange(desc(med)) %>% select(drug)
    drug_levels <- drug_levels$drug
    flt_drug_data$drug <- factor(flt_drug_data$drug,levels=drug_levels)
    p <- ggplot(data=flt_drug_data, aes_string(x="drug", y="AC50", group="sample")) 
    p <- p + geom_point(aes(color=sample), size=3) + theme_bw()
    p + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab('Drug') + ylab('AC50')
  })
  
  #output$drugScreen_QC_plot1 <- renderPlot({
  drugScreen_QC_plot1 <- reactive({
    flt_drug_data <- get_drug_data()
    p <- ggplot(data=flt_drug_data, aes(x=AUC, fill=sample)) + geom_density(alpha=.7)
    p + theme_bw()
  })
  
  #output$drugScreen_QC_plot2 <- renderPlot({
  drugScreen_QC_plot2 <- reactive({
    flt_drug_data <- get_drug_data()
    p <- ggplot(data=flt_drug_data, aes(x=maxResp, fill=sample)) + geom_density(alpha=.7)
    p + theme_bw()
  })
  
  #output$drugScreen_QC_plot3 <- renderPlot({
  drugScreen_QC_plot3 <- reactive({
    validate(need(!all(is.na(data$AC50)), "AC50 data does not exist." ))
    flt_drug_data <- get_drug_data()
    p <- ggplot(data=flt_drug_data, aes(x=AC50, fill=sample)) + geom_density(alpha=.7)
    p + theme_bw()
  })
  
  #output$drugScreen_QC_plot4 <- renderPlot({
  drugScreen_QC_plot4 <- reactive({
    validate(need(!all(is.na(data$IC50)), "IC50 data does not exist." ))
    flt_drug_data <- get_drug_data()
    p <- ggplot(data=flt_drug_data, aes(x=AC50, fill=sample)) + geom_density(alpha=.7)
    p + theme_bw()
  })
  
  QC_plot_list <- reactive({
    plotlist <- list(drugScreen_QC_plot1,drugScreen_QC_plot2)
    if(!all(is.na(data$AC50))){
      plotlist <- c(plotlist,drugScreen_QC_plot3)
    }
    if(!all(is.na(data$IC50))){
      plotlist <- c(plotlist,drugScreen_QC_plot4)
    }
    
    plotlist
  })
  
  output$QC_plots <- renderPlot({
    multiplot(plotlist = QC_plot_list(),cols = 2)
  })
  
  output$drugScreen_dataTable <- renderDataTable(
    get_drug_data()
  )
  
  output$downloadData <- downloadHandler(
    filename = function() { 'demo.csv' },
    content = function(file) {
      write.csv(get_drug_data(), file)
    }
  )
}

# helper function from Abhi
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
