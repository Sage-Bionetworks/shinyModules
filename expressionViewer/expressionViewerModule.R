expressionViewerModuleUI <- function(id){
  ns <- NS(id)
  tagList(
  myHeader <- dashboardHeader(disable=TRUE),
  mySidebar <- dashboardSidebar(disable=TRUE),
  
  myBody <-dashboardBody(
    tags$head(tags$style(HTML('
                                .col-sm-3, {
                                  padding:0 !important;
                                }
                              '))),
    fluidRow(
      column(width = 2,style='padding:0px;',
             
             # Choose sample labels
             box(width=NULL, status='primary', collapsible=TRUE, 
                 collapsed=FALSE, solidHeader=TRUE,
                 title = tagList(shiny::icon("th-list", lib="glyphicon"),
                                 "Label samples"),               
                 uiOutput(ns("anno"))
             ),
             
             # Select genes
             box(width=NULL, status='primary', collapsible=FALSE, 
                 collapsed=FALSE, solidHeader=TRUE,
                 title = tagList(shiny::icon("check", lib="glyphicon"),
                                 "Select genes"),
                 selectInput(ns("select_by"), "Select by", choices = c("gene name","pathway"),selected = "gene name"),
                 uiOutput(ns("genes")),
                 actionButton(ns("refreshGene"), "Refresh")
             ),
             
             #Clustering box
             box(width = NULL, status = "warning", solidHeader=TRUE, 
                 collapsible=TRUE, collapsed=TRUE,
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
      column(width = 10,style='padding:0px;',
             tabBox(width = 12, #solidHeader = TRUE,
                 tabPanel("Heatmap",
                          plotOutput(ns("heatmap"), height = 650)
                          ),
                 tabPanel("PCA",
                          helpText("Always colored by the first selected annotation label."),
                          plotOutput(ns("PCA_plot"))
                          ),
                 tabPanel("Drugs",
                          h5("placeholder")
                          )
             )   
            )    
    )
  )
  )
  dashboardPage(header=myHeader, sidebar=mySidebar, body=myBody,
                skin = "blue")
}

expressionViewerModule <- function(input,output,session,data,tag){
  dataset <- reactive({
    ds <- data
    flog.debug(sprintf("filtered ds dims: %s", dim(ds)), name="server")
    rows_to_keep <- order(apply(exprs(ds),1,var),decreasing=T)[1:50]
    ds_filtered <- ds[rows_to_keep, ]
    
    ds_filtered
  })
  
  metadata <- reactive({
    m_eset <- dataset()
    metaData <- pData(m_eset)
    metaData
  })

  ns <- NS(tag)
  output$anno <- renderUI({
    metaData <- metadata()
    tagList(
      selectInput(ns('annotation_labels'),'Annotate Samples by:',
                  choices=colnames(metaData), selectize=T, multiple=T, selected=colnames(metaData)[1])
    )
  })
  
  output$genes <- renderUI({
    m_eset <- dataset()
    m <- exprs(m_eset)
    geneList <- rownames(m)

    tagList(
      tags$textarea(paste0(c(geneList), collapse="\n"),
                    rows=5, id=ns("selected_genes"), style="width: 100%")
    )
  })
  
  filtered_dataset <- reactive({
    ds <- dataset()
    if(input$refreshGene){
      geneList <- isolate(input$selected_genes)
      geneList <- clean_list(geneList)
      geneList<- intersect(geneList, rownames(fData(ds)))
      ds <- ds[geneList,]
    }
    ds
  })
  
  heatmap_cache <- reactiveValues()
  
  anno_labels <- reactive({
    validate(need(length(input$annotation_labels) > 0, "Please select at least 1 label."))
    validate(need(length(input$annotation_labels) <= 2, "Please select at most 2 labels."))
    input$annotation_labels
  })
  
  #return the heatmap plot
  output$heatmap <- renderPlot({  
    flog.debug("Making heatmap", name='server')
    
    cluster_rows <- input$cluster_rows
    cluster_cols <- input$cluster_cols
    
    m_eset <- filtered_dataset()
    m <- exprs(m_eset)
    m <- data.matrix(m)
    
    validate( need( ncol(m) != 0, "Filtered matrix contains 0 samples.") )
    validate( need( nrow(m) != 0, "Filtered matrix contains 0 genes.") )
    
    metadata <- metadata()
    anno <- anno_labels()
    annotation <- get_heatmapAnnotation(anno, metadata)
    
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
  
  output$PCA_plot <- renderPlot({
    data <- exprs(filtered_dataset())
    m_data <- metadata() 
    anno <- anno_labels()[1]
    pca_res <- prcomp(t(data), center=F, scale=F)
    df <- data.frame(pca_res$x[,c(1:5)])

    percent_variation <- pca_res$sdev^2/sum(pca_res$sdev^2) * 100

    df <- merge(df,m_data,by="row.names")
    rownames(df) <- df$Row.names
    df$Row.names <- NULL

    plotlist <- lapply(1:4, function(i){
      x <- paste0('PC',i)
      y <- paste0('PC',i+1)
      p <- ggplot(data=df, aes_string(x=x,y=y,color=anno))
      p <- p + geom_point() + theme_bw(base_size = 14)
      p <- p + xlab(paste0(x,' - (', round(percent_variation[i],2), '%)' ))  + ylab(paste0(y,' - ( ', round(percent_variation[i+1],2), '%)' ))
    })
    grid_arrange_shared_legend(plotlist, ncol = 2, nrow = 2)
  })
}

grid_arrange_shared_legend <- function(plots, ncol = length(plots), nrow = 1, position = c("bottom", "right")) {
  
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}

clean_list <- function(x) {
  # Split by space, comma or new lines
  x <- unlist(strsplit(x, split=c('[\\s+,\\n+\\r+)]'),perl=T))
  
  # remove the blank entries
  x <- x[!(x == "")]
  
  x
}