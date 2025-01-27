gc()
rm(list = ls())

code.version <- "integrate_BrainMet_datasets"
options(future.globals.maxSize = 10000 * 1024^2)
PROJECT <- "integrated_BrainMet_dataset"

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

##### use the function s8_integration_and_clustering_SeuratV5.selectedGenes.R
source(file.path(path.to.pipeline.src, "processes_src", "s8_integration_and_clustering_SeuratV5.selectedGenes.R"))

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
path.to.14.output <- file.path(path.to.main.output, "14_output")
path.to.16.output <- file.path(path.to.main.output, "16_output")
dir.create(path.to.16.output, showWarnings = FALSE, recursive = TRUE)

#####----------------------------------------------------------------------#####
##### Read data from 12_output
#####----------------------------------------------------------------------#####
num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25
regressOut_mode <- NULL
features_to_regressOut <- NULL
use.sctransform <- TRUE
vars.to.regress <- c("percent.mt")
cluster.resolution <- 0.5

if (file.exists(file.path(path.to.16.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds")) == FALSE){
  # read data from the v0.2 dataset and the GSE193745 dataset
  s.obj.GSE193745 <- readRDS(file.path(path.to.11.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
  s.obj.prev <- readRDS(file.path(path.to.09.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
  
  selected.cells <- c(colnames(s.obj.GSE193745), colnames(s.obj.prev))
  
  ##### get data from v0.3 integration version
  s.obj <- readRDS(file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", "v0.3"), "02_output", "integrated.rds"))
  
  s.obj <- subset(s.obj, cells = selected.cells)
  print(sprintf("Number of cells left after filtering: %s", length(colnames(s.obj))))
  
  ##### Get ribosome genes
  Ribo.gene.patterns <- c("Rpl", "Rps", "RPL", "RPS") 
  all_genes <- row.names(s.obj)
  Ribo.genes <- unlist(lapply(all_genes, function(x){
    if (substr(x, 1, 3) %in% Ribo.gene.patterns){
      return(x)
    }
  }))
  
  num.PCA <- 25
  num.PC.used.in.UMAP <- 25
  num.PC.used.in.Clustering <- 25
  regressOut_mode <- NULL
  features_to_regressOut <- NULL
  use.sctransform <- TRUE
  vars.to.regress <- c("percent.mt")
  cluster.resolution <- 0.5
  
  print("remove samples that have less than 150 cells ...")
  count.cell.in.samples <- table(s.obj$name) %>% data.frame()
  keep.samples <- subset(count.cell.in.samples, count.cell.in.samples$Freq >= 150)$Var1
  
  s.obj <- subset(s.obj, name %in% keep.samples)
  print(sprintf("Number of genes: %s", length(row.names(s.obj))))
  ##### remove ribo genes
  s.obj <- subset(s.obj, features = setdiff(row.names(s.obj), Ribo.genes))
  
  print(sprintf("Number of genes: %s", length(row.names(s.obj))))
  
  DefaultAssay(s.obj) <- "RNA"
  print("Running join layers ...")
  s.obj <- JoinLayers(s.obj)
  
  print("Running integration...")
  
  #> Note that we still save the output to integrated version v0.2
  
  s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = s.obj, 
                                                       save.RDS.s8 = TRUE,
                                                       path.to.output = path.to.16.output,
                                                       use.sctransform = TRUE,
                                                       num.PCA = num.PCA,
                                                       num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                       num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                       cluster.resolution = cluster.resolution,
                                                       vars.to.regress = vars.to.regress,
                                                       remove.genes = NULL)  
  
} else {
  print(sprintf(
    "Data exists at %s",
    file.path(path.to.16.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds")
  ))
  s.obj.integrated <- readRDS(file.path(path.to.16.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds"))
}
