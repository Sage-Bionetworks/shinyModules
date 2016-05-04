#download the DMSO norm data
drug_normViab.synId <- 'syn2773792'
drug_normViab.file <- synGet(drug_normViab.synId)
drug_normViab <- read.delim(drug_normViab.file@filePath, check.names=F, sep="\t", header=T)
drug_normViab$cellLine <- gsub("^ ", "", drug_normViab$cellLine)

#drop unnecassary cols
drop_cols <- c('plate', 'medianDMSO', 'viability')
drug_normViab <- drug_normViab[, !colnames(drug_normViab) %in% drop_cols]

#download the precomputed IC Vals
drug_ICVals.synId <- 'syn2773794'
drug_ICVals.file <- synGet(drug_ICVals.synId)
drug_ICVals <- read.delim(drug_ICVals.file@filePath, check.names=F, sep="\t", header=T)
drug_ICVals <- filter(drug_ICVals, goodNess_of_fit > .70 & hillSlope < 0 )
