# global.R
options(stringsAsFactors = FALSE)
library(shiny)
library(shinydashboard)
library(plyr)
library(synapseClient)

synapseLogin()

source('dataSummarizationModule.R')
source('getData.R')
