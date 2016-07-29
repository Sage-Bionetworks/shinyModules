# example from synergyfinder

## Get combinedData List
combinedData1 <- ReshapeData(mathews_screening_data)
combinedData2 <- combinedData1
combinedData3 <- combinedData1

combinedData <- list(sample1_format1_metric1 = combinedData1, 
                     sample2_format2_metric2 = combinedData2,
                     sample1_format1_metric2 = combinedData3)


# Get sample Informations
sampleNames <- names(combinedData)

sampleInfo <- sapply(sampleNames,function(x){
  parts <- strsplit(x,"_")
  return(unlist(parts))
})
row.names(sampleInfo) <- c("sample","format","metric")
sampleInfo <- as.data.frame(t(sampleInfo))

