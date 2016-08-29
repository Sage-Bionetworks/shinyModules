flog.info('downloading mRNA Exp data from Synapse', name='synapse')
expr_obj <- synGet('syn7124098')

#read in the file
expr <- read.csv(expr_obj@filePath,sep = "\t",stringsAsFactors = FALSE)

flog.info('Reading mRNA metadata from Synapse', name='synapse')
metaData_obj <- synGet("syn7139168")
metaData <- read.csv(metaData_obj@filePath,sep = "\t",stringsAsFactors = FALSE)

## Only keep samples in both
in_common <- intersect(rownames(metaData), colnames(expr))
metaData <- metaData[in_common, ]
expr <- expr[, in_common]
metaData$Synapse.ID <- NULL

features <- data.frame(explicit_rownames=rownames(expr))
rownames(features) <- rownames(expr)

# Scale rows and columns
#expr <- scale(expr)
expr <- t(scale(t(expr)))

eset.data <- ExpressionSet(assayData=as.matrix(expr),
                           phenoData=AnnotatedDataFrame(metaData),
                           featureData=AnnotatedDataFrame(features))


MSIGDB_syn<-synGet("syn2227979")
load(MSIGDB_syn@filePath) #available as MSigDB R object
pathways_list <- c(MSigDB$C2.CP.BIOCARTA, MSigDB$C2.CP.KEGG, MSigDB$C2.CP.REACTOME)
