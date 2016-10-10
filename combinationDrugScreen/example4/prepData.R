library(synapseClient)
synapseLogin()
library("tidyr")


# Synodos NF2 Combination Drug Screen Data
screendirs <- list(`10x10`='syn7210596',`6x6`='syn7210556')

#'Get list of all cells available for combo and measure
#'@param comboScreen
#'@return data frame from synapse table query
getCellForCombo<-function(comboScreen=c('10x10','6x6'))
{
  fileparent=screendirs[[comboScreen]]  
  cells=synapseQuery(paste("select name,id from entity where parentId=='",fileparent,"'",sep=''))
  print(paste('Retrieved',nrow(cells),'cell types for',comboScreen,'screen'))
  return(cells)
}

getData<-function(sample,comboScreen=c('10x10','6x6'),metric=c('CTG','CCG')){
  synid <- cells[match(sample,cells[,1]),2]
  allf <- synapseQuery(paste("select name,id from entity where parentId=='",synid,"'",sep=''))
  
  key_meta <- "metadata"
  key_resp <- "responses"
  # for synodos data only
  if(comboScreen=='10x10'){
    key_meta <- paste("72hours",key_meta,sep="_")
    key_resp <- paste("72hours",key_resp,sep="_")
  }
  meta<-data.frame(fread(synGet(allf[grep(key_meta,allf[,1]),2])@filePath))
  resp<-data.frame(fread(synGet(allf[grep(key_resp,allf[,1]),2])@filePath))
  
  size <- meta$Size[1]
  
  ##first separate out row and column concentrations
  rowCon <- paste('row',c(1:size),sep='_')
  colCon <- paste('col',c(1:size),sep='_')
    
  fullDf<-tidyr::separate(meta,RowConcs,into=rowCon,sep=',')%>%tidyr::separate(ColConcs,into=colCon,sep=',')
  
  rowdf<-fullDf%>%select(match(c("BlockId",rowCon),names(fullDf)))%>%tidyr::gather('Row','RowConc',2:(size+1))%>%arrange(BlockId)
  coldf<-fullDf%>%select(match(c("BlockId",colCon),names(fullDf)))%>%tidyr::gather('Col','ColConc',2:(size+1))%>%arrange(BlockId)
  
  ##now rename rows and columns
  rowdf$Row<-sapply(rowdf$Row,function(x) gsub('row_','',x))
  coldf$Col<-sapply(coldf$Col,function(x) gsub('col_','',x))
  
  metdf<-select(fullDf,match(c('BlockId','Size','RowName','ColName','RowConcUnit','ColConcUnit'),colnames(fullDf)))
  
  ##join rows and columns, sort both the metadata and responses in the same order
  jm<-full_join(metdf,rowdf,by='BlockId')%>%full_join(coldf,by='BlockId')%>%arrange(Col)%>%arrange(Row)%>%arrange(BlockId)
  resp<-resp%>%arrange(Col)%>%arrange(Row)%>%arrange(BlockId)
  
  final.df<-bind_cols(jm,select(resp,c(Value,Replicate)))%>%select(c(BlockId,Size,Row,Col,Value,Replicate,RowName,ColName,RowConc,ColConc,RowConcUnit))
  
  final.df$Row<-as.integer(final.df$Row)
  final.df$Col<-as.integer(final.df$Col)
  final.df$RowConc<-as.numeric(final.df$RowConc)
  final.df$ColConc<-as.numeric(final.df$ColConc)
  final.df<-final.df%>%arrange(Col)%>%arrange(Row)%>%arrange(BlockId)
  
  final.df$sample <- sample
  final.df$assay <- metric
  
  final.df <- final.df[,c("sample","metric","Size","RowName","ColName","RowConc","ColConc","Value","RowConcUnit")]
  colnames(final.df)<-c("sample","assay","numDosagePoints","drug1","drug2","conc1","conc2","response","concUnit")
  return(final.df)
}

getDataList <- function(cells, format = c("10x10","6x6"), metric=c("CTG","CCG")){
  result <- lapply(cells$entity.name, function(x){
    df <- getData(sample = x, comboScreen = format,metric = metric)
    return(df)
  })
  result <- do.call("rbind",result)
  return(result)
}

#10x10 data
format <- "10x10"
cells <- getCellForCombo(format)

metric <- "CTG"
result1 <- getDataList(cells,format,metric)

#6x6 data
format <- "6x6"
cells <- getCellForCombo(format)

metric <- "CTG"
result2 <- getDataList(cells,format,metric)

combinedData <- rbind(result1,result2)

combinedData$drug1 <- editDrugNames(combinedData$drug1)
combinedData$drug2 <- editDrugNames(combinedData$drug2)

editDrugNames <- function(col){
  col <- sub("AMG 232","AMG-232",col)
  col <- sub("AZD2014","AZD-2014",col)
  col <- sub("^Carfilzomib.*","Carfilzomib",col)
  col <- sub("MLN7243","MLN-7243",col)
  return(col)
}

combinedData <- combinedData[combinedData$drug1 != "null" & combinedData$drug2 != "null",]

save(combinedData, file="example4.RData")
