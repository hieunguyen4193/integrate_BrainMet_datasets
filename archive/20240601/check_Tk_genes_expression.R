gc()
rm(list = ls())

code.version <- "20240601"
options(future.globals.maxSize = 10000 * 1024^2)
PROJECT <- "integrated_BrainMet_dataset"

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))
library(ggpubr)
path.to.main.src <- file.path("/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq", code.version)
all.integrated.metadata <- list(v0.1 = read.csv(file.path(path.to.main.src, "samples_to_integrated_20240513.csv")))

integrated.version <- "v0.1"
outdir <- "/media/hieunguyen/HNSD_mini/data/outdir"
main.PROJECT <- "BrainMet_SeuratV5"

path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")
path.to.04.output <- file.path(path.to.main.output, "04_output")
dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)
path.to.save.output <- file.path(path.to.04.output, "check_TK_genes_in_all_samples")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
s.obj <- readRDS(file.path(path.to.03.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
tk.gene.set <- read.csv("/media/hieunguyen/HNHD01/data/UKK/geneSet/TK_genes.csv")$TK_gene
for (gene.id in tk.gene.set){
  p <- VlnPlot(object = s.obj, group.by = "name", features = c(gene.id), pt.size = 0) + theme_pubr()  + theme(legend.position = "right", axis.text.x = element_text(angle = 90))
  ggsave(plot = p, filename = sprintf("gene_%s_in_all_samples.png", gene.id), path = path.to.save.output, device = "png", width = 14, height = 10)
}

