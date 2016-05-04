drugScreenModuleUI <- function(id){
  ns <- NS(id)
  tagList(
           sidebarPanel(
             
             h4('1. Select Cell Lines'), 
             selectInput(ns('cellLines'),NULL, choices = c('ALL', unique(drug_normViab$cellLine)),
                         selectize=T, multiple=T, selected = 'ALL'),
             
             h4('2. Select Drugs'),
             selectInput(ns('selected_drugs'),NULL, choices = c('ALL', unique(drug_normViab$drug)),
                         selectize=T, multiple=T, selected = 'ALL'),
             
             br(), br(),
             h4('Plot settings'),
             selectInput(ns('facet_by'),'Facet by', choices =c('experiment', 'cellLine'),
                         selected = c('experiment'),
                         selectize=T, multiple=T)
           ),
           
           mainPanel(
             tabsetPanel(id="drug_screens", type="tabs",
                         tabPanel("Viability",
                                  checkboxInput(ns("drugViability_heatmap_col_cluster"), label=("cluster drugs (y-axis)"), 
                                                value=TRUE),
                                  hr(),
                                  plotOutput(ns("global_drugViab_heatMap"),height="700px",width="auto",hoverId=NULL),
                                  br(), br(),
                                  helpText("The heatmap above compares the normalized cell viablity(%) across the drug dosages(x-axis) \n
                                           Each y-axis row is a replicate for the selected cellLines and drugs.")
                                  
                                  ),
                         tabPanel("Max Efficacy",
                                  plotOutput(ns("drug_efficacy"),height="700px",width="auto",hoverId=NULL),
                                  helpText("ps: Efficacy is the % of cells that were killed by the given drug")
                         ),
                         tabPanel("ICx",
                                  h4('Select Cell Viability % (ICx)'),
                                  sliderInput(ns('selected_IC_value'),'IC Value', min=10,
                                              max=90, step=10, value=50),
                                  plotOutput(ns("drugScreen_ICx_plot"),height="700px",width="auto",hoverId=NULL)
                         ),
                         tabPanel("Dose Response",
                                  radioButtons(ns("dose_response_plot_splitBy"), label=("split graph by"), 
                                               choices = c('drug', 'cellLine'), selected = 'drug'),
                                  plotOutput(ns("drugResponse_plots"),height="700px",width="auto",hoverId=NULL)
                         )
             )
            )
  )
} 

drugScreenModule <- function(input,output,session,data1,data2){
  get_selected_cellLines <- reactive({
    cellLines <- if('ALL' %in% input$cellLines) unique(data1$cellLine) else input$cellLines
    validate(need(length(cellLines) != 0, "At least one cellLine needs to be selected" ) )
    cellLines
  })
  
  get_selected_drugs <- reactive({
    drugs <- if('ALL' %in% input$selected_drugs) unique(data1$drug) else input$selected_drugs
    validate(need(length(drugs) != 0, "At least one drug needs to be selected"))
    drugs
  })
  
  get_drug_flt_normViab <- reactive({
    flt_Drug_normViab <- filter(data1, drug %in% get_selected_drugs())  
    flt_Drug_normViab <- filter(flt_Drug_normViab, cellLine %in% get_selected_cellLines())  
    flt_Drug_normViab
  })
  
  get_drug_flt_ICVals <- reactive({
    flt_drug_ICVals <- filter(data2, drug %in% get_selected_drugs())  
    flt_drug_ICVals <- filter(flt_drug_ICVals, cellLine %in% get_selected_cellLines())  
    flt_drug_ICVals
  })
  
  output$drugScreen_ICx_plot <- renderPlot({
    flt_drug_ICVals <- get_drug_flt_ICVals()
    ICx <- eval(paste0('IC', input$selected_IC_value))
    #remove NA
    flt_drug_ICVals <- flt_drug_ICVals[! is.na(flt_drug_ICVals[ICx]), ]
    #convert to log10
    flt_drug_ICVals[ICx] <- log10(as.numeric(flt_drug_ICVals[,ICx])) 
    #keep rows where log10 ICx <= 0
    flt_drug_ICVals <- flt_drug_ICVals[flt_drug_ICVals[ICx] <= 0,]
    
    
    drug_levels <- flt_drug_ICVals %>%
      group_by(drug) %>%
      summarise(med=median(IC50, na.rm=T)) %>%
      arrange(desc(med)) %>% select(drug)
    drug_levels <- drug_levels$drug
    flt_drug_ICVals$drug <- factor(flt_drug_ICVals$drug,levels=drug_levels)
    facet_by <- paste(input$facet_by, collapse = ' + ' )
    facet_by <- formula(paste(facet_by, ' ~ .'))
    p <- ggplot(data=flt_drug_ICVals, aes_string(x="drug", y=ICx, group="cellLine")) 
    p <- p + geom_point(aes(color=cellLine), size=3) + theme_bw()
    if( length(input$facet_by) > 0){
      p <- p  + facet_grid(facet_by)    
    }
    p + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab('Drug') + ylab(paste0(ICx, ' (log10 molar conc)'))
  })
  
  output$drug_efficacy <- renderPlot({
    flt_drug_ICVals <- get_drug_flt_ICVals()
    
    drug_levels <- flt_drug_ICVals %>%
      group_by(drug) %>%
      summarise(med=median(maxEfficacy, na.rm=T)) %>%
      arrange(desc(med)) %>% select(drug)
    drug_levels <- drug_levels$drug
    flt_drug_ICVals$drug <- factor(flt_drug_ICVals$drug,levels=drug_levels)
    facet_by <- paste(input$facet_by, collapse = ' + ' )
    facet_by <- formula(paste(facet_by, ' ~ .'))
    p <- ggplot(data=flt_drug_ICVals, aes(x=drug, y=maxEfficacy*100, group=cellLine)) 
    p <- p + geom_point(aes(color=cellLine), size=3) + theme_bw()
    if( length(input$facet_by) > 0){
      p <- p  + facet_grid(facet_by)    
    }
    p <- p + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab('Drug') + ylab('% Efficacy')
    p
  })
  
  
  output$global_drugViab_heatMap <- renderPlot({
    
    validate(need(length(get_selected_cellLines()) != 0, paste0(" Please select cellLine/s")))  
    validate(need(length(input$selected_drugs) !=0 , paste0(" Please select < 5 drugs")))  
    
    x <- get_drug_flt_normViab()
    drugViab_dosages <- dcast(x, experiment+stage+cellLine+drug+replicate ~ conc, value.var="normViability",
                              fun.aggregate = function(x) mean(x))
    
    m <- drugViab_dosages[,-c(1:6)] * 100  #to convert fraction to percentage viability
    rowAnnotation <- drugViab_dosages[,c('drug'),drop=F]
    #m.scaled <- t(scale(t(m)))
    
    #convert colnames to microMolar
    colnames(m) <- as.numeric(colnames(m))*1e+6
    #cluster_rows = if(drugViability_heatmap_col_cluster == TRUE) TRUE else NA
    aheatmap(m,
             scale='none',
             distfun="euclidean",
             Colv=NA,
             Rowv =  if(input$drugViability_heatmap_col_cluster == TRUE) TRUE else NA,
             annRow = rowAnnotation,
             info=TRUE,
             cexRow=0,
             main = 'Cell Viability v/s Drug Dosage(microMolar)',
             sub = 'color signifies cell viability %'
    )
  })
  
  
  output$drugResponse_plots <- renderPlot({
    
    validate(need(length(input$selected_drugs) != 0, paste0(" Please select drug/s (max upto 4)")))  
    validate(need(length(input$selected_drugs) < 5, paste0(" Please select < 5 drugs")))  
    flt_drug_normViab <- get_drug_flt_normViab()
    
    doseResp <- ddply(.data=flt_drug_normViab, .variables = c('drug', 'cellLine', 'experiment'), 
                      .fun = tmp_iterator, .parallel = T)
    
    facet_by <- paste(input$facet_by, collapse = ' + ' )
    facet_by <- formula(paste(facet_by, ' ~ ', input$dose_response_plot_splitBy))
    print(facet_by)
    color_options = c('drug', 'cellLine')
    color_by <- color_options[!color_options %in% input$dose_response_plot_splitBy]
    
    p <- ggplot(data=doseResp, aes_string(x="fittedX", y="fittedY*100", group=color_by))
    p <- p + geom_line(aes_string(color=color_by)) + facet_grid(facet_by) + theme_bw()
    p <- p + geom_hline(aes(yintercept=50), color='grey50', linetype='dashed')
    p <- p + xlab('micromolar conc') + ylab('cell viability %') 
    p <- p + scale_x_continuous(breaks=seq(from=-10,to=-1,by=1),
                                labels = lapply(seq(from=-10,to=-1,by=1), function(x) (10^x)*(1e+6)  ))
    
    p
  })
  
  
}