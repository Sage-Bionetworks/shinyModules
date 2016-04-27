heatmapModuleUI <- function(id){
  ns <- NS(id)
  tagList(
  myHeader <- dashboardHeader(title="", disable=TRUE),
  mySidebar <- dashboardSidebar(disable=TRUE),
  
  myBody <-dashboardBody(
    fluidRow(
      column(width = 3,
             
             # Choose sample labels
             box(width=NULL, status='primary', collapsible=TRUE, 
                 collapsed=FALSE, solidHeader=TRUE,
                 title = tagList(shiny::icon("th-list", lib="glyphicon"),
                                 "Label samples"),               
                 selectInput(ns('heatmap_annotation_labels'),'Annotate Samples by:',
                             choices=colnames(metaData), selectize=T, multiple=T, selected=colnames(metaData)[1])
             ),
             #Clustering box
             box(width = NULL, status = "warning", solidHeader=TRUE, 
                 collapsible=TRUE, collapsed=FALSE,
                 title = tagList(shiny::icon("wrench", lib="glyphicon"), "Change cluster options"),
                 #distance metric
                 selectInput(ns("clustering_distance"), "Distance Calculation",
                             choices=c("correlation", "euclidean", "maximum", 
                                       "manhattan", "canberra", "binary", "minkowski"),
                             selectize=T, multiple=F, selected="euclidean"),
                 # set the clustering method
                 selectInput(ns("clustering_method"), "Clustering Method",
                             choices=c("ward", "single", "complete", "average", 
                                       "mcquitty", "median", "centroid"),
                             selectize=T, multiple=F, selected="average"),
                 checkboxInput(ns('cluster_cols'), 'Cluster the columns', value = TRUE),
                 checkboxInput(ns('cluster_rows'), 'Cluster the rows', value = TRUE)
             )
      ),
      column(width = 9,
             box(width = NULL, solidHeader = TRUE,
                 plotOutput(ns("heatmap"), height = 650))
            )    
    )
  )
  )
  dashboardPage(header=myHeader, sidebar=mySidebar, body=myBody,
                skin = "blue")
}

heatmapModule <- function(input,output,session,data){
  filtered_dataset <- reactive({
    ds <- data
    flog.debug(sprintf("filtered ds dims: %s", dim(ds)), name="server")
    rows_to_keep <- order(apply(exprs(ds),1,var),decreasing=T)[1:50]
    ds_filtered <- ds[rows_to_keep, ]
    
    ds_filtered
  })
  
  heatmap_cache <- reactiveValues()
  
  #return the heatmap plot
  output$heatmap <- renderPlot({  
    flog.debug("Making heatmap", name='server')
    
    cluster_rows <- input$cluster_rows
    cluster_cols <- input$cluster_cols
    
    m_eset <- filtered_dataset()
    m <- exprs(m_eset)
    m <- data.matrix(m)
    
    validate( need( ncol(m) != 0, "Filtered matrix contains 0 Samples.") )
    validate( need( nrow(m) != 0, "Filtered matrix contains 0 features.") )
    validate( need(nrow(m) < 10000, "Filtered matrix contains > 10000 genes.") )
    
    filtered_metadata <- pData(m_eset)
    annotation <- get_heatmapAnnotation(input$heatmap_annotation_labels, filtered_metadata)
    
    fontsize_row <- ifelse(nrow(m) > 100, 0, 8)
    fontsize_col <- ifelse(ncol(m) > 50, 0, 8)    
    
    heatmap.color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
    
    
    heatmap_cache$heatmap <- sbHeatMap(m,annotation,
                                        clustering_distance_rows = input$clustering_distance,
                                        clustering_distance_cols = input$clustering_distance,
                                        fontsize_col=fontsize_col, 
                                        fontsize_row=fontsize_row,
                                        scale=F,
                                        color=heatmap.color,
                                        #breaks=heatmap.breaks,
                                        clustering_method = input$clustering_method,
                                        explicit_rownames = fData(m_eset)$explicit_rownames,
                                        cluster_rows=cluster_rows, cluster_cols=cluster_cols,
                                        drawColD=FALSE)
    
  })
}