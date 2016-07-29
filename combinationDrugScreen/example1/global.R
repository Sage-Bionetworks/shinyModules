#global.R
options(stringsAsFactors = FALSE)
library(shiny)
library(shinydashboard)
library("ggplot2")
library("data.table")
library("grid")
library("futile.logger")
library("synergyfinder")

flog.threshold(DEBUG, name='server')
flog.threshold(DEBUG, name='ui')
flog.threshold(DEBUG, name='global')

flog.debug("Starting App...", name="server")

flog.debug("Loading module...", name="server")
source("../combinationDrugScreenModule.R")

flog.debug("Loading data...", name="server")

source("getData.R")
