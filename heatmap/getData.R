
flog.info('Demo: NF2 primary meningioma normalized mRNA Exp data from Synapse', name='synapse')
expr_obj <- synGet('syn5907232')

#read in the file
expr <- read.csv(expr_obj@filePath,sep = "\t",stringsAsFactors = FALSE)


flog.info('Demo: Reading mRNA metadata from Synapse', name='synapse')
metaData_obj <- synGet("syn5987186")
metaData <- read.csv(metaData_obj@filePath,sep = "\t",stringsAsFactors = FALSE)
row.names(metaData) <- metaData$SampleID

## Only keep samples in both
in_common <- intersect(rownames(metaData), colnames(expr))
metaData <- metaData[in_common, ]
expr <- expr[, in_common]

features <- data.frame(explicit_rownames=rownames(expr))
rownames(features) <- rownames(expr)

# Scale rows and columns
expr <- scale(expr)
expr <- t(scale(t(expr)))

eset.data <- ExpressionSet(assayData=as.matrix(expr),
                           phenoData=AnnotatedDataFrame(metaData),
                           featureData=AnnotatedDataFrame(features))
