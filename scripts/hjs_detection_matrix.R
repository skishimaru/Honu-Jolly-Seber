# Creates detection matrix for the multi-state jolly-seber model for honu 
# Created by: Shelbie Ishimaru
# Created on: 2024-09-09
# Last edited: 2024-09-09
################################################################################
# Initialization ---------------------------------------------------------------
library(tidyverse) #for data manipulation
library(here) # for creating unbreakable file paths

# Read-in Data -----------------------------------------------------------------
nesting23 <- read.csv(here("Data", "Nesting23.csv")) #read-in raw nesting data 

# Clean Data ----------------------------------------------------------------
nesting23 <- nesting23[nesting23$Area== "Tern",] #removes all East Island nesting instances

