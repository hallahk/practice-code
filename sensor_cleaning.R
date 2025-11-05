
knitr::opts_chunk$set(echo = TRUE)


# Purple Air data cleaning and calibrating process: co-location process

# Clear environment
rm(list = ls(all = TRUE))
if (!is.null(sessionInfo()$otherPkgs)) {
  res <- suppressWarnings(
    lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
           detach, character.only=TRUE, unload=TRUE, force=TRUE))
}


# Load packages, installing if needed
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
pacman::p_load(
  riem,
  openair,
  scales,
  RColorBrewer,
  table1,
  plyr,
  dplyr,
  tidyr,
  tidyverse,
  ggplot2,
  MetBrewer,
  stringr,
  data.table,
  ggpubr,
  zoo,
  tidyquant,
  here
)


## Checking whether PurpleAirs are performing well:


# Functions for reading in co-location datasets and getting sensor names
read_and_name_csv <- function(file_path) {
  # Extract folder and file name without extension
  folder_name <- dirname(file_path) %>% basename()
  file_name <- basename(file_path) %>% tools::file_path_sans_ext()
  
  # Create the name for the data frame
  df_name <- paste(folder_name, file_name, sep = "_")
  
  # Read the CSV file
  df <- read.csv(file_path)
  
  # Assign the data frame to the global environment with the constructed name
  assign(df_name, df, envir = .GlobalEnv)
}

#Get sensor name function from the folder name
get_sensor_name <- function(df_name) {
  sensor_name <- str_sub(df_name, start = 18, end = -14)
  sensor_name
}

