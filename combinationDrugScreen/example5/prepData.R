library(synapseClient)
synapseLogin()
library("tidyr")
library("data.table")
library("dplyr")

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

getDataHelper <- function(sample,assay,allf){
  
  meta<-data.frame(fread(synGet(allf[grep("metadata",allf[,1]),2])@filePath))
  resp<-data.frame(fread(synGet(allf[grep("responses",allf[,1]),2])@filePath))
  calc<-data.frame(fread(synGet(allf[grep("calc",allf[,1]),2])@filePath)) 
  
  meta_calc <- merge(meta,calc[,c("BlockId","DBSumNeg")],by="BlockId")
  
  size <- meta_calc$Size[1]
  
  ##first separate out row and column concentrations
  rowCon <- paste('row',c(1:size),sep='_')
  colCon <- paste('col',c(1:size),sep='_')
  
  fullDf<-tidyr::separate(meta_calc,RowConcs,into=rowCon,sep=',')%>%tidyr::separate(ColConcs,into=colCon,sep=',')
  
  rowdf<-fullDf%>%select(match(c("BlockId",rowCon),names(fullDf)))%>%tidyr::gather('Row','RowConc',2:(size+1))%>%arrange(BlockId)
  coldf<-fullDf%>%select(match(c("BlockId",colCon),names(fullDf)))%>%tidyr::gather('Col','ColConc',2:(size+1))%>%arrange(BlockId)
  
  ##now rename rows and columns
  rowdf$Row<-sapply(rowdf$Row,function(x) gsub('row_','',x))
  coldf$Col<-sapply(coldf$Col,function(x) gsub('col_','',x))
  
  metdf<-select(fullDf,match(c('BlockId','Size','RowName','ColName','RowConcUnit','ColConcUnit','DBSumNeg'),colnames(fullDf)))
  
  ##join rows and columns, sort both the metadata and responses in the same order
  jm<-full_join(metdf,rowdf,by='BlockId')%>%full_join(coldf,by='BlockId')%>%arrange(Col)%>%arrange(Row)%>%arrange(BlockId)
  resp<-resp%>%arrange(Col)%>%arrange(Row)%>%arrange(BlockId)
  
  final.df<-bind_cols(jm,select(resp,c(Value,Replicate)))%>%select(c(BlockId,Size,Row,Col,Value,Replicate,RowName,ColName,RowConc,ColConc,RowConcUnit,DBSumNeg))
  
  final.df$RowConc<-as.numeric(final.df$RowConc)
  final.df$ColConc<-as.numeric(final.df$ColConc)
  final.df<-final.df%>%arrange(Col)%>%arrange(Row)%>%arrange(BlockId)
  
  final.df$sample <- sample
  final.df$assay <- assay
  
  final.df <- final.df[,c("sample","assay","Size","RowName","ColName","RowConc","ColConc","Value","RowConcUnit","DBSumNeg")]
  colnames(final.df)<-c("sample","assay","numDosagePoints","drug1","drug2","conc1","conc2","response","concUnit","DBSumNeg")
  return(final.df)
}

getData<-function(sample,comboScreen=c('10x10','6x6'),assay=c('CTG','CCG')){
  synid <- cells[match(sample,cells[,1]),2]
  allf <- synapseQuery(paste("select name,id from entity where parentId=='",synid,"'",sep=''))
  if(comboScreen == '10x10'){
    allf<-allf[grep(assay,allf[,1]),]
  }
  if(length(allf[,1]) == 6){
    set1f <- allf[grep("set1",allf[,1]),]
    df1 <- getDataHelper(sample,assay,set1f)
    set2f <- allf[grep("set2",allf[,1]),]
    df2 <- getDataHelper(sample,assay,set2f)
    df <- rbind(df1,df2)
  }else{
    df <- getDataHelper(sample,assay,allf)
  }
  return(df)
}

getDataList <- function(cells, format = c("10x10","6x6"), assay=c("CTG","CCG")){
  result <- lapply(cells$entity.name, function(x){
    df <- getData(sample = x, comboScreen = format,assay = assay)
    return(df)
  })
  result <- do.call("rbind",result)
  return(result)
}

#10x10 data
format <- "10x10"
cells <- getCellForCombo(format)

assay <- "CTG"
result1 <- getDataList(cells,format,assay)

assay <- "CCG"
result2 <- getDataList(cells,format,assay)

#6x6 data
format <- "6x6"
cells <- getCellForCombo(format)

assay <- "CTG"
result3 <- getDataList(cells,format,assay)

combinedData <- do.call(rbind,list(result1,result2,result3))

editDrugNames <- function(col){
  col <- sub("\\s+$","",col)
  col <- sub("AMG-232","AMG 232",col)
  col <- sub("AZD-2014","AZD2014",col)
  col <- sub("^BKM-120$","BKM120 (NVP-BKM120, Buparlisib)",col)
  col <- sub("^Cabozantinib$","Cabozantinib (XL184, BMS-907351)",col)
  col <- sub("^Carfilzomib$","Carfilzomib (PR-171)",col)
  col <- sub("^Crizotinib$","Crizotinib (PF-02341066)",col)
  col <- sub("^Danusertib$","Danusertib (PHA-739358)",col)
  col <- sub("^Doxorubicin$","Doxorubicin (Adriamycin)",col)
  col <- sub(".+JQ1","JQ1",col)
  col <- sub("^LY2801653$","LY2801653 dihydrochloride",col)
  col <- sub("MK-2206 2HCl","MK-2206",col)
  col <- sub("MLN-7243","MLN7243",col)
  col <- sub("Topotecan HCl","Topotecan hydrochloride",col)
  col <- sub("^PD-0332991$","Palbociclib (PD-0332991) HCl",col)
  return(col)
}


combinedData$drug1 <- editDrugNames(combinedData$drug1)
combinedData$drug2 <- editDrugNames(combinedData$drug2)


combinedData <- combinedData[combinedData$drug1 != "null" & combinedData$drug2 != "null",]
combinedData <- combinedData[combinedData$drug1 != "N/A ",]

save(combinedData, file="example5.RData")
