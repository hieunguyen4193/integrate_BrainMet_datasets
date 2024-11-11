gc()
rm(list = ls())

library(Seurat)
input.data <- Read10X_h5("/media/hieunguyen/HNHD01/data/UKK/BrainMet/raw_downloaded_data/GSE174401/GSE174401_filtered_feature_bc_matrix.h5") 
saveRDS(input.data, file.path("/media/hieunguyen/HNHD01/data/UKK/BrainMet/preprocess_raw_data/GSE174401/GSE174401/GSE174401.rds"))
