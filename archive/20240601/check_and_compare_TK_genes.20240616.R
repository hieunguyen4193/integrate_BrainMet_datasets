gc()
rm(list = ls())

code.version <- "20240601"
options(future.globals.maxSize = 10000 * 1024^2)
PROJECT <- "integrated_BrainMet_dataset"

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
source(file.path(path.to.pipeline.src, "processes_src", "s8_integration_and_clustering_SeuratV5.R"))
source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))
library(Matrix)
path.to.main.src <- file.path("/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq", code.version)
all.integrated.metadata <- list(v0.1 = read.csv(file.path(path.to.main.src, "samples_to_integrated_20240513.csv")))

integrated.version <- "v0.1"
outdir <- "/media/hieunguyen/HNSD_mini/data/outdir"
main.PROJECT <- "BrainMet_SeuratV5"

path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")
dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)
s.obj.integrated <- readRDS(file.path(path.to.03.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))

sample.metadata <- read.csv("/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq/20240601/samples_to_integrated_20240513.csv")
meta.data <- s.obj.integrated@meta.data %>% rownames_to_column("barcode")
meta.data$group <- unlist(lapply(meta.data$name, function(x){
  return(subset(sample.metadata, sample.metadata$Sample == x)$Group)}))
meta.data <- meta.data %>% column_to_rownames("barcode")
meta.data <- meta.data[row.names(s.obj.integrated@meta.data), ]
s.obj.integrated <- AddMetaData(object = s.obj.integrated, metadata = meta.data$group, col.name = "group")

tk.gene.set <- read.csv("/media/hieunguyen/HNHD01/data/UKK/geneSet/TK_genes.csv")$TK_gene

# FeaturePlot(object = s.obj.integrated, reduction = "harmony_UMAP", features = c("CSF1R", "CD3D", "CD3E", "LCK", "FYN", "ZAP70"), label = TRUE, repel = TRUE, ncol = 3 )
# 
# FeaturePlot(object = s.obj.integrated, reduction = "harmony_UMAP", features = c("PTPRC"), label = TRUE, repel = TRUE)
# 
# VlnPlot(object = s.obj.integrated, features = c("CSF1R", "CD3D", "CD3E", "LCK", "FYN", "ZAP70"))
# 
# to_vec( for (item in row.names(s.obj.integrated)) if(grepl("PTPRC", item) == TRUE) item)

meta.data <- s.obj.integrated@meta.data %>% rownames_to_column("barcode")

sample.cells <- list()
for (dataset.name in unique(meta.data$name)){
  sample.cells[[dataset.name]] <- subset(meta.data, meta.data$name == dataset.name)$barcode
}

groups <- list()
for (input.group in unique(meta.data$group)){
  groups[[input.group]] <- unique(subset(meta.data, meta.data$group == input.group)$name)
}

checkdf <- data.frame(Gene = tk.gene.set)
all.count.mat <- GetAssayData(object = s.obj.integrated, slot = "data", assay = "SCT")[tk.gene.set, ]

get_pct <- function(dataset.name, gene.id){
  count.mat <- all.count.mat[gene.id, sample.cells[[dataset.name]]]
  return(length(count.mat[count.mat != 0])/length(count.mat))  
}

for (dataset.name in unique(s.obj.integrated$name)){
  print(sprintf("working on %s", dataset.name))
  checkdf[[dataset.name]] <- unlist(lapply(checkdf$Gene, function(x){
    return(get_pct(dataset.name, x))
  }))
}

library(heatmaply)
checkdf <- checkdf %>% column_to_rownames("Gene")
p.all <- ggheatmap(checkdf, scale = "none",
               row_dend_left = FALSE, row_text_angle = 0,
               column_text_angle = 90,
               scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                 low = "blue",
                 high = "red"))


p.brainmet <- ggheatmap(checkdf[, groups$BrainMet], scale = "none",
                   row_dend_left = FALSE, row_text_angle = 0,
                   column_text_angle = 90,
                   scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                     low = "blue",
                     high = "red"))


p.control<- ggheatmap(checkdf[, groups$control], scale = "none",
                   row_dend_left = FALSE, row_text_angle = 0,
                   column_text_angle = 90,
                   scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                     low = "blue",
                     high = "red"))
