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

colnames(integrated.metadata) <- c("PROJECT", "dataset.name", "sample.id")
integrated.metadata <- integrated.metadata %>% rowwise() %>%
  mutate(input.name = sprintf("%s_%s_%s", PROJECT, dataset.name, sample.id))

outdir <- "/media/hieunguyen/HNSD_mini/data/outdir"
main.PROJECT <- "BrainMet_SeuratV5"

path.to.main.output <- file.path(outdir, main.PROJECT, sprintf("integrated_%s", integrated.version))
path.to.02.output <- file.path(path.to.main.output, "02_output")
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)

PROJECT <- "integrated_BrainMet_dataset"
input.outdir <- "/media/hieunguyen/HNHD01/outdir"

data.list <- list()

for (i in seq(1, nrow(integrated.metadata))){
  input.name <- integrated.metadata[i,]$input.name
  dataset.name <- integrated.metadata[i, ]$dataset.name
  sample.id <- integrated.metadata[i, ]$sample.id

  path.to.main.input <- file.path(input.outdir, input.name, sprintf("%s_%s", dataset.name, sample.id), "s1_output")
  input.file <- Sys.glob(file.path(path.to.main.input, sprintf("*%s*", dataset.name)))
  if (length(input.file) != 1){
    print(sprintf("Error at %s", sample.id))
  } else {
    print(sprintf("adding sample %s", input.name))
    data.list[[sample.id]] <- readRDS(input.file[[1]])
    meta.data <- data.list[[sample.id]]@meta.data %>% rownames_to_column("barcode")
    meta.data$sample.id <- sample.id
    meta.data <- meta.data %>% column_to_rownames("barcode")
    meta.data <- meta.data[row.names(data.list[[sample.id]]@meta.data),]
    data.list[[sample.id]] <- AddMetaData(object = data.list[[sample.id]], metadata = meta.data$sample.id, col.name = "name")
    print(sprintf("dataset size: %s", paste(dim(data.list[[sample.id]]), collapse = ", ")))
  }
}

print("Finished merging...")
s.obj <- merge(data.list[[1]], data.list[2:length(data.list)])
print("Saving...")
saveRDS(s.obj, file.path(path.to.02.output, "integrated.rds"))

print("Finished gemerating merged input for integration process.")