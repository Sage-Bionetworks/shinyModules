library("Biobase")
flog.info('Demo: Generating Random Gene Expression data for 100 genes and 50 samples')

#read in the file
ngenes <- 100
nsamples <- 50
expr <- matrix(rnorm(ngenes*nsamples), nrow=ngenes)
geneNames <- paste0('gene-', 1:nrow(expr) )
rownames(expr) <- geneNames

sampleNames <- paste0('sample-', 1:ncol(expr))
colnames(expr) <- sampleNames

features <- data.frame(explicit_rownames=rownames(expr))
rownames(features) <- rownames(expr)

#introducing bias 
expr[,1:30] <- expr[,1:30] + 10    # high disease
expr[,31:50] <- expr[,31:50] + 4  # females
metaData <- data.frame('age'= sample(10:80,ncol(expr)),
                       'gender'=c(rep('male',30),rep('female',20)),
                       'disease'=c(rep('high',30), rep('low',20)))
rownames(metaData) <- sampleNames


# Scale rows and columns
m <- t(scale(t(expr)))

eset.data <- ExpressionSet(assayData=as.matrix(m),
                           phenoData=AnnotatedDataFrame(metaData),
                           featureData=AnnotatedDataFrame(features))

MSIGDB_syn<-synGet("syn2227979")
load(MSIGDB_syn@filePath) #available as MSigDB R object
pathways_list <- c(MSigDB$C2.CP.BIOCARTA, MSigDB$C2.CP.KEGG, MSigDB$C2.CP.REACTOME)
