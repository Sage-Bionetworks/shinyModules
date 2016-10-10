#global.R
options(stringsAsFactors = FALSE)
library(shiny)
library(shinydashboard)
library("ggplot2")
library("data.table")
library("gridExtra")
library("futile.logger")
library("dplyr")
library("reshape")
library(plyr)
library(doMC)
registerDoMC(4)
library(nplr)
library(stringr)

flog.threshold(DEBUG, name='server')
flog.threshold(DEBUG, name='ui')
flog.threshold(DEBUG, name='global')
flog.threshold(INFO, name='synapse')

flog.debug("Starting App...", name="server")

flog.debug("Loading module...", name="server")
source("../combinationDrugScreenModule.R")

flog.debug("Loading data...", name="server")

load("example4.RData")
