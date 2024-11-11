gc()
rm(list = ls())

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))
library(Matrix)
path.to.main.src <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq"
all.integrated.metadata <- list(v0.1 = read.csv(file.path(path.to.main.src, "samples_to_integrated_20240513.csv")))

integrated.version <- "v0.1"

#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
integrated.metadata <- all.integrated.metadata[[integrated.version]]

colnames(integrated.metadata) <- c("PROJECT", "dataset.name", "SampleID")
integrated.metadata <- integrated.metadata %>% rowwise() %>%
  mutate(input.name = sprintf("%s_%s_%s", PROJECT, dataset.name, SampleID))

outdir <- "/media/hieunguyen/HNSD_mini/data/outdir"
main.PROJECT <- "BrainMet_SeuratV5"

path.to.main.output <- file.path(outdir, main.PROJECT, sprintf("integrated_%s", integrated.version))
path.to.02.output <- file.path(path.to.main.output, "02_output")
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)

PROJECT <- "integrated_BrainMet_dataset"
input.outdir <- "/media/hieunguyen/HNHD01/outdir"

#####----------------------------------------------------------------------#####
##### Step 1: Generate matrix file .rds
#####----------------------------------------------------------------------#####
# for (i in seq(1, nrow(integrated.metadata))){
#   input.name <- integrated.metadata[i,]$input.name
#   dataset.name <- integrated.metadata[i, ]$dataset.name
#   sampleid <- integrated.metadata[i, ]$SampleID
#   
#   path.to.main.input <- file.path(input.outdir, input.name, sprintf("%s_%s", dataset.name, sampleid), "data_analysis/01_output")
#   input.file <- Sys.glob(file.path(path.to.main.input, sprintf("*%s*", dataset.name)))
#   if (length(input.file) != 1){
#     print(sampleid)
#   } else {
#     print(sprintf("adding sample %s", input.name))
#     tmp.s.obj <- readRDS(input.file)
#     tmp.mat <- GetAssayData(tmp.s.obj, slot = "counts", assay = "RNA")
#     saveRDS(tmp.mat, file.path(path.to.02.output, sprintf("mat_sample_%s.rds", input.name)))
#     rm(tmp.mat)
#   }  
#   gc()
# }

#####----------------------------------------------------------------------#####
##### Step 2: Merge matrix file into one file
#####----------------------------------------------------------------------#####
# dir.create(file.path(path.to.02.output, "merge_mat"), showWarnings = FALSE, recursive = TRUE)
# 
# files <- Sys.glob(file.path(path.to.02.output, "*.rds"))
# splitted.files <- split(files, seq(1,5))
# 
# for (part.idx in names(splitted.files)){
#   if (file.exists(file.path(path.to.02.output, "merge_mat", sprintf("mat.part%s.rds", part.idx))) == FALSE){
#     all.files <- splitted.files[[part.idx]]
#     mat <- readRDS(all.files[[1]])  
#     finished.samples <- c(all.files[[1]])
#     for (file in all.files[2:length(all.files)]){
#       print(sprintf("WORKING ON SAMPLE %s", basename(file)))
#       tmpdf <- readRDS(file)
#       print(sprintf("Shape: %s", paste(dim(tmpdf), collapse = ", ")))
#       merge.mat <- merge(mat, tmpdf, by="row.names", all = TRUE) %>% column_to_rownames("Row.names")
#       merge.mat[is.na(merge.mat)] <- 0
#       mat <- Matrix(as.matrix(merge.mat), sparse = TRUE)
#       rm(merge.mat)
#       gc()
#       print(sprintf("Shape of merged matrix: %s", paste(dim(mat), collapse = ", ")))
#       finished.samples <- c(finished.samples, file)
#     }
#     saveRDS(mat, file.path(path.to.02.output, "merge_mat", sprintf("mat.part%s.rds", part.idx)))
#     writexl::write_xlsx(data.frame(status = finished.samples), file.path(path.to.02.output, "merge_mat", sprintf("finished_samples.part%s.xlsx",part.idx)))
#   }
#   rm(mat)
#   gc()
# }

#####----------------------------------------------------------------------#####
##### Step 3: Merge .part matrix file into one file
#####----------------------------------------------------------------------#####

# files <- Sys.glob(file.path(path.to.02.output, "merge_mat", "mat.part*.rds"))
# mat <- list()
# mat.colnames <- list()
# data.list <- list()
# 
# for (i in seq(1, length(files))){
#   mat[[sprintf("s%s", i)]] <- readRDS(files[[i]])  
#   mat.colnames[[i]] <- colnames(mat[[sprintf("s%s", i)]])
#   data.list[[i]] <- CreateSeuratObject(counts = mat[[sprintf("s%s", i)]])
# }
# 
# s.obj <- merge(data.list[[1]], data.list[2:length(data.list)])
# s.obj <- JoinLayers(s.obj)
# saveRDS(GetAssayData(object = s.obj, slot = "counts", assay = "RNA"), file.path(path.to.02.output, "integrated_mat.rds"))

#####----------------------------------------------------------------------#####
##### RE-DO
#####----------------------------------------------------------------------#####
# data.list <- list()
# all.single.mat <- Sys.glob(file.path(path.to.02.output, "single_mat", "*.rds"))
# for (file in all.single.mat){
#   sample.id <- basename(file)
#   sample.id <- str_replace(sample.id, "mat_sample_BrainMet_SeuratV5_", "")
#   sample.id <- str_replace(sample.id, ".rds", "")
#   print(sprintf("adding sample %s", sample.id))
#   data.list[[sample.id]] <- CreateSeuratObject(counts = readRDS(file))
#   meta.data <- data.list[[sample.id]]@meta.data %>% rownames_to_column("barcode")
#   meta.data$sample.id <- sample.id
#   meta.data <- meta.data %>% column_to_rownames("barcode")
#   meta.data <- meta.data[row.names(data.list[[sample.id]]@meta.data),]
#   data.list[[sample.id]] <- AddMetaData(object = data.list[[sample.id]], metadata = meta.data$sample.id, col.name = "sample.id")
# }
# 
# s.obj.integrated <- merge(data.list[[1]], data.list[2:length(all.single.mat)])
# saveRDS(s.obj.integrated, file.path(path.to.02.output, "integrated_dataset", "integrated.rds"))
