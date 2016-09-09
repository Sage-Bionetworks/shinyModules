#global.R
options(stringsAsFactors = FALSE)
library(shiny)
library(shinydashboard)
library(synapseClient)
library("ggplot2")
library("data.table")
library("grid")
library("futile.logger")
library("synergyfinder")
library("tidyr")
library("dplyr")

flog.threshold(DEBUG, name='server')
flog.threshold(DEBUG, name='ui')
flog.threshold(DEBUG, name='global')
flog.threshold(INFO, name='synapse')

synapseLogin()

flog.debug("Starting App...", name="server")

flog.debug("Loading module...", name="server")
source("../combinationDrugScreenModule2.R")

flog.debug("Loading data...", name="server")

#source("getData.R")
load("example3_norm.RData")
