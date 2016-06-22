options(stringsAsFactors = F)
options(warn=-1)

library("devtools")
library("synapseClient")
library("shiny")
library("shinydashboard")
library("dplyr")
library("tidyr")
library("memoise")
library("futile.logger")
library("Biobase")
library("data.table")
library("ggplot2")

# Set up logging
flog.threshold(DEBUG, name='server')
flog.threshold(DEBUG, name='ui')
flog.threshold(DEBUG, name='global')
flog.threshold(INFO, name='synapse')

#login to synapse
synapseLogin()

flog.debug("Starting App", name="server")

source("../../lib/sbHeatmap.R")

#get module
source("../expressionViewerModule.R")

#get data
source("getData.R")