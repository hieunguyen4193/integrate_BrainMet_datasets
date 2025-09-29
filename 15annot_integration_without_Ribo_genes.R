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
path.to.15.output <- file.path(path.to.main.output, "15_output")
dir.create(path.to.15.output, showWarnings = FALSE, recursive = TRUE)

path.to.save.output <- file.path(file.path(path.to.15.output, sprintf("ribo_thres_%s", "INTEGRATE_WITHOUT_RIBO_GENES")) )

s.obj <- readRDS(file.path(path.to.save.output, "integrated_BrainMet_dataset.moreClusterRes.output.s8.rds"))
annot.metadata <- read.csv("/media/hieunguyen/GSHD_HN01/outdir/BrainMet_SeuratV5/integrate_BrainMet_datasets/integrated_v0.2/integration_evaluation/metadata_barcode_celltype.csv")
annot.metadata <- subset(annot.metadata, select = c(barcode, celltype))

path.to.project.src <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq/integrate_BrainMet_datasets"
sample.metadata <- read.csv(file.path(path.to.project.src, "samples_to_integrated_20240725.csv"))
meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
meta.data <- merge(meta.data, sample.metadata, by.x = "name", by.y = "Sample")
meta.data <- meta.data %>% column_to_rownames("barcode")
meta.data <- meta.data[row.names(s.obj@meta.data), ]

s.obj <- AddMetaData(object = s.obj, col.name = "label", metadata = meta.data$Group)

barcode.list <- list()
for (input.label in unique(s.obj$label)){
  barcode.list[[input.label]] <- colnames(subset(s.obj, label == input.label))
}
meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
meta.data <- merge(meta.data, annot.metadata, by.x = "barcode", by.y = "barcode", all.x = TRUE)
meta.data <- meta.data %>% column_to_rownames("barcode")
meta.data <- meta.data[row.names(s.obj@meta.data), ]
s.obj <- AddMetaData(object = s.obj, metadata = meta.data$celltype, col.name = "celltype")

##### subset and plot
# selected.celltype <- "TAMs"
# selected.label <- "Control Glioma"

dir.create(file.path(path.to.save.output, "UMAP_celltype_label"), showWarnings = FALSE, recursive = TRUE)
for (selected.label in unique(s.obj$label)){
  subset.s.obj <- subset(s.obj, cells = barcode.list[[selected.label]])
  for (selected.celltype in unique(subset.s.obj$celltype)){
    p <- DimPlot(object = subset(subset.s.obj, celltype == selected.celltype) , 
                 reduction = "harmony_UMAP", 
                 label = TRUE, 
                 label.box = TRUE, 
                 repel = TRUE, 
                 group.by = "celltype")    
    ggsave(plot = p, filename = sprintf("%s.%s.pdf", selected.celltype, selected.label), 
           path = file.path(path.to.save.output, "UMAP_celltype_label"), 
           device = "pdf", 
           width = 14, 
           height = 10, 
           dpi = 300)
  }
}

