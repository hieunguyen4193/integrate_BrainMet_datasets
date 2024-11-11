gc()
rm(list = ls())

library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)

maindir <- "/media/hieunguyen/HNHD01/data/UKK/BrainMet"
dataset.name <- "GSE131907"

path.to.rds.file <- file.path(maindir, "raw_downloaded_data", dataset.name, "GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
path.to.save.output <- file.path(maindir, "preprocess_raw_data")
input.data <- readRDS(path.to.rds.file)

barcodedf <- data.frame(barcode = colnames(input.data))
barcodedf <- barcodedf %>% rowwise() %>%
  mutate(sample_name = paste(str_split(barcode, "_")[[1]][2:3] , collapse = "_"))

for (sample.id in unique(barcodedf$sample_name)){
  print(sprintf("working on sample %s", sample.id))
  dir.create(file.path(path.to.save.output, sprintf("%s_%s", dataset.name, sample.id)), showWarnings = FALSE, recursive = TRUE)
  selected.barcodes <- subset(barcodedf, barcodedf$sample_name == sample.id)$barcode
  output.rds <- input.data[, selected.barcodes]
  saveRDS(output.rds, file.path(path.to.save.output, sprintf("%s_%s", dataset.name, sample.id), sprintf("%s_%s", dataset.name, sample.id)))
}
