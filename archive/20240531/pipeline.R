gc()
rm(list = ls())

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
path.to.main.src <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq"

# source(file.path(path.to.main.src, "02_merging_samples.R"))
source(file.path(path.to.main.src, "03_manually_run_integration.R"))