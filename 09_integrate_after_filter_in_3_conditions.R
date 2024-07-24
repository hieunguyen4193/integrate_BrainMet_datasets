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

outdir <- "/media/hieunguyen/HNSD01/outdir"
main.PROJECT <- "BrainMet_SeuratV5"

sample.metadata <- all.integrated.metadata[[integrated.version]]
path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")
path.to.09.output <- file.path(path.to.main.output, "09_output")
dir.create(path.to.09.output, showWarnings = FALSE, recursive = TRUE)

running.cases <- list(
  BrainMet = list(
    cluster.resolution = 0.5,
    remove.cluster = c(15)),
  `Control Epilepsy` = list(
    cluster.resolution = 0.75,
    remove.cluster = c(5, 7)
  ),
  `Control Glioma` = list(
    cluster.resolution = 0.5,
    remove.cluster = c(10, 11)
  )
)


if (file.exists(file.path(path.to.09.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT))) == FALSE){
  
  selected.cells <- c()
  for ( i in seq(1, length(running.cases))){
    input.condition <- names(running.cases)[[i]]
    remove.clusters <- running.cases[[input.condition]][["remove.cluster"]]
    cluster.resolution <- running.cases[[input.condition]][["cluster.resolution"]]
    path.to.08.output <- file.path(path.to.main.output, sprintf("08_output_%s_%s", input.condition, cluster.resolution))
    tmp <- readRDS(file.path(path.to.08.output, "s8_output/integrated_BrainMet_dataset.output.s8.rds"))
    selected.cells <- c(selected.cells, colnames(tmp))
  }
  
  s.obj <- readRDS(file.path(path.to.02.output, "integrated.rds"))
  s.obj <- subset(s.obj, cells = selected.cells)
  
  num.PCA <- 25
  num.PC.used.in.UMAP <- 25
  num.PC.used.in.Clustering <- 25
  regressOut_mode <- NULL
  features_to_regressOut <- NULL
  use.sctransform <- TRUE
  vars.to.regress <- c("percent.mt")
  cluster.resolution <- 0.5
  
  DefaultAssay(s.obj) <- "RNA"
  s.obj <- JoinLayers(s.obj)
  
  s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = s.obj, 
                                                       save.RDS.s8 = TRUE,
                                                       path.to.output = path.to.09.output,
                                                       use.sctransform = TRUE,
                                                       num.PCA = num.PCA,
                                                       num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                       num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                       cluster.resolution = cluster.resolution,
                                                       vars.to.regress = vars.to.regress)    
} else {
  print("data exists, reading in ...")
  s.obj.integrated <- readRDS(file.path(path.to.09.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
}

