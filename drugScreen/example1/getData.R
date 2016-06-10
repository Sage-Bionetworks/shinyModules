
#NCATS 
synodosData <- synGet("syn5816995")
synodosData <- read.xls(synodosData@filePath)
names(synodosData) <- c("protocolName","drugID","sample","sampleType","nf2Status",
                        "AC50","curveClass","maxResp","logAC50","drug","AUC",
                        "AUCfit","target")

select_col <- c("sample","AC50","maxResp","curveClass","drugID","drug","target","AUC")
summarizedData <- synodosData[,select_col]

summarizedData$IC50 <- NA
