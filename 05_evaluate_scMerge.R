gc()
rm(list = ls())

code.version <- "integrate_BrainMet_datasets"
options(future.globals.maxSize = 10000 * 1024^2)
PROJECT <- "integrated_BrainMet_dataset"

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

path.to.main.src <- file.path("/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq", code.version)
all.integrated.metadata <- list(v0.1 = read.csv(file.path(path.to.main.src, "samples_to_integrated_20240513.csv")))

integrated.version <- "v0.2"
outdir <- "/media/hieunguyen/HNSD01/outdir"
main.PROJECT <- "BrainMet_SeuratV5"

cluster.resolution <- 0.5

path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))
path.to.06.output <- file.path(path.to.main.output, sprintf("06_output_%s", cluster.resolution))

path.to.sobj <- file.path(path.to.06.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT))
s.obj <- readRDS(path.to.sobj)
meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
  mutate(batch = name) %>%
  column_to_rownames("barcode")

meta.data <- meta.data[row.names(s.obj@meta.data),]
s.obj <- AddMetaData(object = s.obj, metadata = meta.data$batch, col.name = "batch")
# if ("scMerge" %in% installed.packages() == FALSE){
#   BiocManager::install("scMerge", update = FALSE)
# }
library(scMerge)
library(scater)

##### selecting negative control, stably expressed genes
data("segList_ensemblGeneID", package = "scMerge") 
data("segList", package = "scMerge")
seg.genes <- segList$human$human_scSEG

s.obj <- RunUMAP(object = s.obj, features = intersect(seg.genes, row.names(s.obj)), reduction.name = "UMAP_SEG", slot = "data", assay = "SCT")
DimPlot(object = s.obj, reduction = "UMAP_SEG", label = TRUE, label.box = TRUE, group.by = "name")
DimPlot(object = s.obj, reduction = "harmony_UMAP", label = TRUE, label.box = TRUE, group.by = "name")
FeaturePlot(object = s.obj, features = sample(seg.genes, 9), ncol = 3, reduction = "harmony_UMAP")

# sce.obj <- as.SingleCellExperiment(s.obj)
# 
# num.batch <- length(unique(s.obj$batch))
# scMerge_unsupervised <- scMerge(
#   sce_combine = sce.obj,
#   kmeansK = rep(3, num.batch),
#   ctl = seg.genes,
#   assay_name = "scMerge_unsupervised")
# 
# scMerge_unsupervised = runPCA(scMerge_unsupervised, exprs_values = "scMerge_unsupervised")
# scater::plotPCA(
#   scMerge_unsupervised, 
#   colour_by = "batch")
