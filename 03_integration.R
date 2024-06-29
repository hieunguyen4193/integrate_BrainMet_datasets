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

path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")
dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)

s.obj <- readRDS(file.path(path.to.02.output, "integrated.rds"))

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

count.cell.in.samples <- table(s.obj$name) %>% data.frame()
keep.samples <- subset(count.cell.in.samples, count.cell.in.samples$Freq >= 150)$Var1

if (file.exists(file.path(path.to.03.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT))) == FALSE){
  s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = subset(s.obj, name %in% keep.samples), 
                                                       save.RDS.s8 = TRUE,
                                                       path.to.output = path.to.03.output,
                                                       use.sctransform = TRUE,
                                                       num.PCA = num.PCA,
                                                       num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                       num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                       cluster.resolution = cluster.resolution,
                                                       vars.to.regress = vars.to.regress)    
} else {
  print("data exists, reading in ...")
  s.obj.integrated <- readRDS(file.path(path.to.03.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
}

min.pct <- 0.5
Idents(s.obj.integrated) <- "cca.cluster.0.5"
if (file.exists(file.path(path.to.03.output, sprintf("DE_cluster_marker_genes_min_pct_%s.rds", min.pct))) == FALSE){
  print("Start preparation for find all markers ...")
  DefaultAssay(s.obj.integrated) <- "SCT"
  s.obj.integrated <- PrepSCTFindMarkers(s.obj.integrated)
  print("start running FindAllMarkers ...")
  cluster.markers <- FindAllMarkers(object = s.obj.integrated, assay = "SCT", test.use = "wilcox", slot = "data", min.pct = min.pct)
  cluster.markers <- subset(cluster.markers, cluster.markers$p_val_adj < 0.05 & cluster.markers$avg_log2FC > 0)
  saveRDS(cluster.markers, file.path(path.to.03.output, sprintf("DE_cluster_marker_genes_min_pct_%s.rds", min.pct)))
} else {
  cluster.markers <- readRDS(file.path(path.to.03.output, sprintf("DE_cluster_marker_genes_min_pct_%s.rds", min.pct)))
}
