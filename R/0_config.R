# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
# load packages
library('ctmle')
library('data.table')
library('dplyr')
library('foreach')
library('tmle3')
library('sl3')
library('dplyr')
library('stringr')
library('glmnet')
library('here')
library('tidyverse')
library('forcats')
library('origami')
library('hal9001')
library('tictoc')

# subetting helper function
"%w/o%" <- function(x, y) x[!x %in% y]