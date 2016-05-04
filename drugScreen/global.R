#global.R
options(stringsAsFactors = FALSE)
library(shiny)
library(shinyIncubator)
library(synapseClient)
library('rCharts')
library("RCurl")
library("reshape2")
library("scales")
library("gdata")
library("plyr")
library("dplyr")
library("devtools")
library("ggplot2")
library("data.table")
library("doMC")
library("NMF")
registerDoMC(4)

synapseLogin()

source("global_DrugScreens.R")

source("drugScreenModule.R")

source("getData.R")

