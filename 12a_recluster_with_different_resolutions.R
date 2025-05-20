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
dir.create(path.to.12a.output, showWarnings = FALSE, recursive = TRUE)

print("Reading main seurat object data ...")
s.obj <- readRDS(file.path(path.to.12.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds"))

selected.reduction <- "harmony"
num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25
regressOut_mode <- NULL
features_to_regressOut <- NULL
use.sctransform <- TRUE
vars.to.regress <- c("percent.mt")

if (file.exists(file.path(path.to.12a.output, "integrated_BrainMet_dataset.moreClusterRes.output.s8.rds")) == FALSE){
  for (cluster.resolution in c(0.1, 0.2, 0.25, 0.3, 0.4)){
    print("-----------------------------------------------------------")
    print(sprintf("Working on clustering with resolution %s", cluster.resolution))
    print("-----------------------------------------------------------")
    s.obj <- FindNeighbors(s.obj, dims = 1:num.PC.used.in.Clustering, reduction = selected.reduction)
    s.obj <- FindClusters(s.obj, resolution = cluster.resolution, cluster.name = sprintf("%s.cluster.%s", selected.reduction, cluster.resolution))  
  }
  saveRDS(s.obj, file.path(path.to.12a.output, "integrated_BrainMet_dataset.moreClusterRes.output.s8.rds")) 
} else {
  print("reading in saved data ...")
  s.obj <- readRDS(file.path(path.to.12a.output, "integrated_BrainMet_dataset.moreClusterRes.output.s8.rds"))
}

