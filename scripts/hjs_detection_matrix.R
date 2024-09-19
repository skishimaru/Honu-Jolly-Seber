# Creates detection matrix for the multi-state jolly-seber model for honu 
# Created by: Shelbie Ishimaru
# Created on: 2024-09-09
# Last edited: 2024-09-16
################################################################################
# Initialization ---------------------------------------------------------------
library(tidyverse) #for data manipulation
library(tidyr) #for data frame manipulation
library(here) #for creating unbreakable file paths

# Read-in Data -----------------------------------------------------------------
nesting23 <- read.csv(here("Data", "Nesting23.csv")) #read-in raw nesting data 

# Clean Data -------------------------------------------------------------------
nesting23 <- nesting23[nesting23$Area== "Tern",] #removes all East Island nesting instances

nesting23$dt <- mdy(nesting23$EVENT_DATE) #stands for month, day, year, make it a true date
nesting23$jday <- yday(nesting23$dt) #creates column for Julian Day
nesting23$week <- format(nesting23$dt, "%V") #creates column with week as a decimal number

nesting23$TURTLEID <- as.factor(nesting23$TURTLEID) #factor turtle ID col for matrix building
nesting23$jday <- as.factor(nesting23$jday) #factor jday col for matrix building
nesting23$week <- as.factor(nesting23$week) #factor week col for matrix building

# Nightly Detection Matrix -----------------------------------------------------
night_detect <- nesting23 %>% group_by(TURTLEID, jday, .drop= FALSE) %>% #select the only columns needed for analysis
  summarise("count"= n(), "detect"= as.numeric(any(count>= 1))) %>% #calculate the number of times each turtle was seen a night and use counts for detection col
  select(TURTLEID, jday, detect) %>% spread(jday, detect) #remove count col and edit df so it's in detection matrix format

write.table(night_detect, here("outputs", "nightly_detection.txt")) #export nightly detection matrix to csv

# Weekly Detection Matrix ------------------------------------------------------
weekly_detect <- nesting23 %>% group_by(TURTLEID, week, .drop= FALSE) %>% #select the only columns needed for analysis
  summarise("count"= n(), "detect"= as.numeric(any(count>= 1))) %>% #calculate the number of times each turtle was seen a week and use counts for detection col
  select(TURTLEID, week, detect) %>% spread(week, detect) #remove count col and edit df so it's in detection matrix format
  
write.table(weekly_detect, here("outputs", "weekly_detection.txt")) #export weekly detection matrix to csv
