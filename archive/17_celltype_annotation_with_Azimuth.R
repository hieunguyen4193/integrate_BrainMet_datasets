gc()
rm(list = ls())

code.version <- "integrate_BrainMet_datasets"
options(future.globals.maxSize = 10000 * 1024^2)
PROJECT <- "integrated_BrainMet_dataset"

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
source(file.path(path.to.pipeline.src, "processes_src", "s8_integration_and_clustering_SeuratV5.R"))
source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))
library(Matrix)
path.to.main.src <- file.path("/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq", code.version)
all.integrated.metadata <- list(v0.1 = read.csv(file.path(path.to.main.src, "samples_to_integrated_20240513.csv")),
                                v0.2 = read.csv(file.path(path.to.main.src, "samples_to_integrated_20240620.csv")))

integrated.version <- "v0.2"

# outdir <- "/media/hieunguyen/HNSD01/outdir"
outdir <- "/media/hieunguyen/GSHD_HN01/outdir"
main.PROJECT <- "BrainMet_SeuratV5"

path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))
path.to.09.output <- file.path(path.to.main.output, "09_output")
path.to.10.output <- file.path(path.to.main.output, "10_output")
path.to.11.output <- file.path(path.to.main.output, "11_output")
path.to.12.output <- file.path(path.to.main.output, "12_output")
path.to.12a.output <- file.path(path.to.main.output, "12a_output")

use.output <- "12a_output"

if (use.output == "12_output"){
  s.obj <- readRDS(file.path(path.to.12.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds"))  
} else if (use.output = "12a_output"){
  
}

path.to.17.output <- file.path(path.to.main.output, "17_output", "symphony", sprintf("from_%s", use.output))
dir.create(path.to.17.output, showWarnings = FALSE, recursive = TRUE)

# sanity check the data
DimPlot(object = s.obj, reduction = "harmony_UMAP", label = TRUE, label.box = TRUE, group.by = "harmony.cluster.0.5")

# if ("Azimuth" %in% installed.packages() == FALSE){
  BiocManager::install(c("EnsDb.Hsapiens.v86c"), update = FALSE)
  install.packages('/media/hieunguyen/GSHD_HN01/storage/offline_pkgs/BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz', type = "source", repos = NULL)
  remotes::install_github('satijalab/azimuth', ref = 'master', upgrade = FALSE)
# }