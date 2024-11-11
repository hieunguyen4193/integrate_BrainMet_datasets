gc()
rm(list = ls())

code.version <- "20240601"

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))
library(Matrix)
path.to.main.src <- file.path("/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq", code.version)
all.integrated.metadata <- list(v0.1 = read.csv(file.path(path.to.main.src, "samples_to_integrated_20240513.csv")),
                                v0.2 = read.csv(file.path(path.to.main.src, "samples_to_integrated_20240620.csv")))

# integrated.version <- "v0.1"
integrated.version <- "v0.2"

#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
integrated.metadata <- all.integrated.metadata[[integrated.version]]

colnames(integrated.metadata) <- c("PROJECT", "dataset.name", "sample.id", "group")
integrated.metadata <- integrated.metadata %>% rowwise() %>%
  mutate(input.name = sprintf("%s_%s_%s", PROJECT, dataset.name, sample.id))

outdir <- "/media/hieunguyen/HNSD01/outdir"
main.PROJECT <- "BrainMet_SeuratV5"

path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))
path.to.02.output <- file.path(path.to.main.output, "02_output")
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)

PROJECT <- "integrated_BrainMet_dataset"
input.outdir <- "/media/hieunguyen/HNHD01/outdir"

for (i in seq(1, nrow(integrated.metadata))){
  input.name <- integrated.metadata[i,]$input.name
  dataset.name <- integrated.metadata[i, ]$dataset.name
  sample.id <- integrated.metadata[i, ]$sample.id
  
  path.to.main.input <- file.path(input.outdir, input.name, sprintf("%s_%s", dataset.name, sample.id), "s3_output")
  input.file <- Sys.glob(file.path(path.to.main.input, sprintf("*%s*", dataset.name)))
  if (length(input.file) != 1){
    print(sprintf("Error at %s", sample.id))
  }
}