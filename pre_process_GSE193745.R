gc()
rm(list = ls())

library(Seurat)
library(dplyr)
library(tidyverse)
library(stringr)

maindir <-"/media/hieunguyen/HNHD01/outdir/BrainMet_SeuratV5_GSE193745_GSE193745/GSE193745_GSE193745/s8a_output"
savedir <- "/media/hieunguyen/HNHD01/outdir"

s.obj <- readRDS(file.path(maindir, "BrainMet_SeuratV5_GSE193745_GSE193745.output.s8a.rds"))
download.meta.data <- read.csv(file.path(maindir, "GSE193745_Mets2022.TIL.metadata.txt"), sep = "\t") %>%
  rownames_to_column("barcode")

##### NOTE
# split the big dataset into sample-wise dataset
# save to the same folder with the same structure, which can be loaded easily 
# by the merge object function in 02 and integrate in 03.

# One of the sample, LB3935T, was excluded from the combined counts matrix and 
# subsequent integrated analysis, due to low cell number (<=100 cells).
# link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193745
#####

meta.data <- s.obj@meta.data %>% rownames_to_column("rowname") %>%
  mutate(barcode = str_replace(rowname, "GSE193745_GSE193745_", ""))

meta.data <- merge(meta.data, subset(download.meta.data, select = c(barcode, ID, CancerType)), by.x = "barcode", by.y = "barcode")

meta.data <- meta.data %>% rowwise() %>%
  mutate(SampleID = ifelse(grepl("_", ID) == TRUE, str_split(ID, "_")[[1]][[1]], ID)) %>%
  column_to_rownames("rowname")

meta.data <- meta.data[row.names(s.obj@meta.data), ]
s.obj <- AddMetaData(object = s.obj, metadata = meta.data$CancerType, col.name = "PrimaryTumor")
s.obj <- AddMetaData(object = s.obj, metadata = meta.data$SampleID, col.name = "SampleID")

dataset.metadata <- data.frame()
for (sample.id in unique(s.obj$SampleID)){
  dir.create(file.path(savedir, sprintf("BrainMet_SeuratV5_GSE193745_%s", sample.id), sprintf("GSE193745_%s", sample.id), "s8_output"), showWarnings = FALSE, recursive = TRUE)
  print(sprintf("saving %s ...", sample.id))
  tmp.s.obj <- subset(s.obj, SampleID == sample.id)
  saveRDS(tmp.s.obj, file.path(savedir, sprintf("BrainMet_SeuratV5_GSE193745_%s", sample.id), sprintf("GSE193745_%s", sample.id), "s8_output", sprintf("BrainMet_SeuratV5_GSE193745_%s.output.s8a.rds", sample.id)))
  
  tmp.dataset.metadata <- data.frame(PROJECT = "BrainMet_SeuratV5",
                                     Dataset = "GSE193745",
                                     Sample = sample.id,
                                     Group = "BrainMet")
  dataset.metadata <- rbind(dataset.metadata, tmp.dataset.metadata)
}
write.csv(dataset.metadata, file.path(maindir, "GSE193745_sample_metadata.csv"))






