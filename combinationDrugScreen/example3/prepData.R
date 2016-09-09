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


#'Format synergy data to be used by the synergyfinder R package
#'There is a bug in the synergy calculators so try to use https://github.com/sgosline/synergyfinder
#'@param cell cell of interest
#'@param comboScreen either 6x6 or 10x10
#'@param measure either CTG or CCG (for 10x10 only)
#'@return list of data to be analyzed by synergyfinder
#'
getSynFinderData<-function(cell,comboScreen=c('10x10','6x6'),measure=c('CTG','CCG')){
#   if(measure=='CCG'&&comboScreen=='6x6'){
#     print('CCG data not available for 6x6, evaluating with CTG')
#     measure='CTG'
#   }
  
  
  #  cell.dat<-lapply(cells[,1],function(x){
  
  synid=cells[match(cell,cells[,1]),2]
  allf<-synapseQuery(paste("select name,id from entity where parentId=='",synid,"'",sep=''))
  #allf<-allf[grep(measure,allf[,1]),]
  key_meta <- "metadata"
  key_resp <- "responses"
  if(comboScreen=='10x10'){
    key_meta <- paste("72hours",key_meta,sep="_")
    key_resp <- paste("72hours",key_resp,sep="_")
  }
  meta<-data.frame(fread(synGet(allf[grep(key_meta,allf[,1]),2])@filePath))
  resp<-data.frame(fread(synGet(allf[grep(key_resp,allf[,1]),2])@filePath))
  
  ##first separate out row and column concentrations
  if(comboScreen=='6x6'){
    rowCon=paste('row',c(1,2,3,4,5,6),sep='_')
    colCon=paste('col',c(1,2,3,4,5,6),sep='_')
  }else{
    rowCon=paste('row',c(1,2,3,4,5,6,7,8,9,10),sep='_')
    colCon=paste('col',c(1,2,3,4,5,6,7,8,9,10),sep='_')
  }
  fullDf<-tidyr::separate(meta,RowConcs,into=rowCon,sep=',')%>%tidyr::separate(ColConcs,into=colCon,sep=',')
  
  if(comboScreen=='6x6'){
    rowdf<-fullDf%>%select(match(c("BlockId",rowCon),names(fullDf)))%>%tidyr::gather('Row','RowConc',2:7)%>%arrange(BlockId)
    coldf<-fullDf%>%select(match(c("BlockId",colCon),names(fullDf)))%>%tidyr::gather('Col','ColConc',2:7)%>%arrange(BlockId)
  } else{
    rowdf<-fullDf%>%select(match(c("BlockId",rowCon),names(fullDf)))%>%tidyr::gather('Row','RowConc',2:11)%>%arrange(BlockId)
    coldf<-fullDf%>%select(match(c("BlockId",colCon),names(fullDf)))%>%tidyr::gather('Col','ColConc',2:11)%>%arrange(BlockId)
    
  }   
  ##now rename rows and columns
  rowdf$Row<-sapply(rowdf$Row,function(x) gsub('row_','',x))
  coldf$Col<-sapply(coldf$Col,function(x) gsub('col_','',x))
  
  metdf<-select(fullDf,match(c('BlockId','RowName','ColName','RowConcUnit','ColConcUnit'),colnames(fullDf)))
  
  ##join rows and columns, sort both the metadata and responses in the same order
  jm<-full_join(metdf,rowdf,by='BlockId')%>%full_join(coldf,by='BlockId')%>%arrange(Col)%>%arrange(Row)%>%arrange(BlockId)
  resp<-resp%>%arrange(Col)%>%arrange(Row)%>%arrange(BlockId)
  
  final.df<-bind_cols(jm,select(resp,c(Value,Replicate)))%>%select(c(BlockId,Row,Col,Value,Replicate,RowName,ColName,RowConc,ColConc,RowConcUnit))
  colnames(final.df)<-c('BlockID','Row','Col','Response','Replicate','DrugRow','DrugCol','ConcRow','ConcCol','ConcUnit')
  final.df$ConcRow<-as.numeric(final.df$ConcRow)
  final.df$ConcCol<-as.numeric(final.df$ConcCol)
  ##return rowtrug, row target, col drug, col target, and some synergy score....
  return(final.df)
}

renameCol <- function(x){
  paste(x,format,metric,sep = "_")
}

getDataList <- function(cells, format = c("10x10","6x6"), metric=c("CTG","CCG")){
  result <- lapply(cells$entity.name, function(x){
    df <- getSynFinderData(cell = x, comboScreen = format,measure = metric)
    result <- ReshapeData(df)
    return(result)
  })
  names(result) <- vapply(cells$entity.name,renameCol,"")
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

combinedData <- c(result1,result2)

# Get sample Informations
sampleNames <- names(combinedData)

sampleInfo <- sapply(sampleNames,function(x){
  parts <- strsplit(x,"_")
  return(unlist(parts))
})
row.names(sampleInfo) <- c("sample","format","metric")
sampleInfo <- as.data.frame(t(sampleInfo))

normData <- lapply(combinedData,function(x){
  temp <- x$drug.pairs
  temp$drug.row <- editDrugNames(temp$drug.row)
  temp$drug.col <- editDrugNames(temp$drug.col)
  return(list("dose.response.mats" = x$dose.response.mats,
              "drug.pairs" = temp))
})

editDrugNames <- function(col){
  col <- sub("AMG 232","AMG-232",col)
  col <- sub("AZD2014","AZD-2014",col)
  col <- sub("^Carfilzomib.*","Carfilzomib",col)
  col <- sub("MLN7243","MLN-7243",col)
  return(col)
}

save(combinedData,sampleInfo, file="example3.RData")
save(normData,sampleInfo,file="example3_norm.RData")
