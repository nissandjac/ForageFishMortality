## Load the data ##
load('data/sprat_datras.Rdata') # Sprat data from Datras

# Load sprat data from the SMS model and the official assessment
Sprat.obs <- read.csv('data/Sprat_sms.csv')
sprat.ass <- read.table('data/sprat_assessment.txt', header = T)
Sprat.M <- read.csv('data/spratNS_mortality.csv')
Sprat.weight <- read.csv('data/spratNS_weight.csv')

# Load required files 
source('gradient.R')
source('baseparameters_trimmed.R')
library(ggplot2)
library(dplyr)
library(scales)


