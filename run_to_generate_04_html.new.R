gc()
rm(list = ls())

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

output.version <- "20240708"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq/integrate_BrainMet_datasets"
path.to.rmd <- file.path(path.to.project.src, "04_downstream_analysis_separate_conditions.Rmd")

all.integrated.metadata <- list(v0.1 = read.csv(file.path(path.to.project.src, "samples_to_integrated_20240513.csv")),
                                v0.2 = read.csv(file.path(path.to.project.src, "samples_to_integrated_20240620.csv")))
integrated.version <- "v0.2"
outdir <- "/media/hieunguyen/HNSD01/outdir"
main.PROJECT <- "BrainMet_SeuratV5"

sample.metadata <- all.integrated.metadata[[integrated.version]]

main.PROJECT <- "BrainMet_SeuratV5"
code.version <- "integrate_BrainMet_datasets"
PROJECT <- "integrated_BrainMet_dataset"

path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))

