library(synapseClient)
synapseLogin()
library(plyr)
library(dplyr)
library(data.table)
library(tidyr)



raw_drugresponse_matrix <- fread(synGet("syn5613652")@filePath, sep=",")
metadata <- fread(synGet("syn5613651")@filePath, sep=",")

View(raw_drugresponse_matrix)


block1 <- raw_drugresponse_matrix %>% filter(BlockId ==  1)
View(block1)
block1 <- block1 %>% spread(Col, Value,)
block1$BlockId <- NULL
block1$Replicate <- NULL
rownames(block1) <- block1$Row
block1$Row <- NULL
image(as.matrix(block1))

image(block1)

View(metadata)


View(block1)

View(block1)

View(block1)
colnames(raw_drugresponse_matrix)
View(raw_drugresponse_matrix)
unique(raw_drugresponse_matrix$BlockId)


stable(raw_drugresponse_matrix$Col)
table(raw_drugresponse_matrix$Row)
table(raw_drugresponse_matrix$Replicate)
