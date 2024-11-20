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

outdir <- "/media/hieunguyen/GSHD_HN01/outdir"
main.PROJECT <- "BrainMet_SeuratV5"

path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))
path.to.03.output <- file.path(path.to.main.output, "03_output")
path.to.09.output <- file.path(path.to.main.output, "09_output")
path.to.10.output <- file.path(path.to.main.output, "10_output")
path.to.11.output <- file.path(path.to.main.output, "11_output")
path.to.12.output <- file.path(path.to.main.output, "12_output")
path.to.13.output <- file.path(path.to.main.output, "13_output")
dir.create(path.to.13.output, showWarnings = FALSE, recursive = TRUE)

s.obj.integrated <- readRDS(file.path(path.to.12.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds"))

##### Update 19.11.2024
gene.list <- c("CD3E", "CD4", "CD8A", "NCAM1", "ITGAX", "CECAM8")
DefaultAssay(s.obj.integrated) <- "SCT"
feature.plot <- FeaturePlot(object = s.obj.integrated, reduction = "harmony_UMAP", features = gene.list, label = TRUE, pt.size = 1) &
  scale_color_gradient(low = "lightgray", high = "#FF0000", na.value = "lightgray")
# scale_color_distiller(palette = "RdBu")
violin.plot <- VlnPlot(object = s.obj.integrated, features = gene.list, pt.size = 0, group.by = "harmony.cluster.0.5")
s.obj.integrated <- AddModuleScore(object = s.obj.integrated, features = list(gene.list), name = "check_gene_list_", ctrl = 50)

feature.plot.module <- FeaturePlot(object = s.obj.integrated, reduction = "harmony_UMAP", features = "check_gene_list_1", label = TRUE, pt.size = 1) &
  scale_color_gradient(low = "lightgray", high = "#FF0000", na.value = "lightgray")
violin.plot.module <- VlnPlot(object = s.obj.integrated, features = "check_gene_list_1", pt.size = 0, group.by = "harmony.cluster.0.5")
