#getData.R

projectDf <- data.frame(projectName = c('Drug Screening of pNF Cell Lines','pNF Cell Line Characterization Studies'), 
                        annotationTable = c('syn7506024','syn7805075'))

keyList <- c('assay','dataType','dataSubtype','fileFormat','species','diagnosis','tumorType',
             'cellType','isCellLine','organ','tissue','consortium')