gc()
rm(list = ls())
# 
library(comprehenr)
library(dplyr)
library(tidyverse)
library(stringr)
# path.to.main.input <- "/media/hieunguyen/HNHD01/data/UKK/BrainMet/input_data"
# 
# all.samples <- Sys.glob(file.path(path.to.main.input, "*", "*barcodes.tsv*"))
# path.to.all.input.data <- "/media/hieunguyen/HNHD01/data/UKK/BrainMet/input_data/all_input_data"
# 
# for (file in all.samples){
#   print(sprintf("working on sample %s", file))
#   sample.name <- str_split(basename(file), "_")[[1]][[1]]
#   sample.id <- sprintf("%s_%s", basename(dirname(file)), sample.name) 
#   dir.create(file.path(path.to.all.input.data, sample.id), showWarnings = FALSE, recursive = TRUE)
#   path.to.save.sample <- file.path(path.to.all.input.data, sample.id)
#   system(sprintf("mv %s %s", file.path(dirname(file), sprintf("%s*", sample.name)), path.to.save.sample))
#   
#   convert.filename <- list(
#     barcodes = "barcodes.tsv.gz",
#     features = "features.tsv.gz",
#     matrix = "matrix.mtx.gz"
#   )
#   for (input.pattern in c("barcodes", "features", "matrix")){
#     all.files <- Sys.glob(file.path(path.to.save.sample, sprintf("*%s.*", input.pattern)))
#     if (length(all.files) == 1){
#       system(sprintf("mv %s %s", all.files[[1]], file.path(dirname(all.files[[1]]), convert.filename[[input.pattern]])))
#     }  
#   }
#   
# }

path.to.main.input <- "/media/hieunguyen/HNHD01/data/UKK/BrainMet/input_data/all_input_data"
all.samples <- Sys.glob(file.path(path.to.main.input, "*"))

for (file in all.samples){
  dir.create(file.path(file, basename(file)), showWarnings = FALSE, recursive = TRUE)
  system(sprintf("mv %s/*.* %s", file, file.path(file, basename(file)) ))
}