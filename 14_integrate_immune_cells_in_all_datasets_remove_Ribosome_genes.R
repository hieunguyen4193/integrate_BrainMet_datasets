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
path.to.14.output <- file.path(path.to.main.output, "14_output")
dir.create(path.to.14.output, showWarnings = FALSE, recursive = TRUE)

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

if (file.exists(file.path(path.to.12.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds")) == FALSE){
  # read data from the v0.2 dataset and the GSE193745 dataset
  s.obj.GSE193745 <- readRDS(file.path(path.to.11.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
  s.obj.prev <- readRDS(file.path(path.to.09.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
  
  selected.cells <- c(colnames(s.obj.GSE193745), colnames(s.obj.prev))
  
  ##### get data from v0.3 integration version
  s.obj <- readRDS(file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", "v0.3"), "02_output", "integrated.rds"))
  
  s.obj <- subset(s.obj, cells = selected.cells)
  print(sprintf("Number of cells left after filtering: %s", length(colnames(s.obj))))
  
  print("remove samples that have less than 150 cells ...")
  count.cell.in.samples <- table(s.obj$name) %>% data.frame()
  keep.samples <- subset(count.cell.in.samples, count.cell.in.samples$Freq >= 150)$Var1
  
  s.obj <- subset(s.obj, name %in% keep.samples)
  
  DefaultAssay(s.obj) <- "RNA"
  print("Running join layers ...")
  s.obj <- JoinLayers(s.obj)
  
  print("Running integration...")
  
  #> Note that we still save the output to integrated version v0.2
  
  s.obj.integrated12 <- s8.integration.and.clustering_V5(s.obj = s.obj, 
                                                       save.RDS.s8 = TRUE,
                                                       path.to.output = path.to.12.output,
                                                       use.sctransform = TRUE,
                                                       num.PCA = num.PCA,
                                                       num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                       num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                       cluster.resolution = cluster.resolution,
                                                       vars.to.regressfeature = vars.to.regress)  
  
} else {
  print(sprintf(
    "Data exists at %s",
    file.path(path.to.12.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds")
  ))
  s.obj.integrated12 <- readRDS(file.path(path.to.12.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds"))
}

#####----------------------------------------------------------------------#####
##### Remove cells by RIBOSOME percentage
#####----------------------------------------------------------------------#####


DefaultAssay(s.obj.integrated12) <- "SCT"
Idents(s.obj.integrated12) <- "harmony.cluster.0.5"

ribo.umap <- FeaturePlot(object = s.obj.integrated12, reduction = "harmony_UMAP", 
                         features = c("percent.ribo"), label = TRUE)
ribo.violin <- VlnPlot(object = s.obj.integrated12, features = c("percent.ribo"), group.by = "harmony.cluster.0.5")

cluster.markers <- readRDS(file.path(path.to.12.output, "s8_output", "DE_cluster_marker_genes_min_pct_0.5_reduction_method_harmony.rds"))

all.ribo.thres <- c(33, 42)

for (ribo.thres in all.ribo.thres){
  ##### rerun the integration after filtering cells
  if (file.exists(file.path(path.to.14.output, 
                            sprintf("ribo_thres_%s", ribo.thres), 
                            "s8_output",
                            "integrated_BrainMet_dataset.output.s8.rds")) == FALSE){
    s.obj <- subset(s.obj.integrated12, percent.ribo < ribo.thres)
    DefaultAssay(s.obj) <- "RNA"
    print("Running join layers ...")
    s.obj <- JoinLayers(s.obj)
    
    print("Running integration...")
    dir.create(file.path(path.to.14.output, sprintf("ribo_thres_%s", ribo.thres)), showWarnings = FALSE, recursive = TRUE)
    s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = s.obj, 
                                                         save.RDS.s8 = TRUE,
                                                         path.to.output = file.path(path.to.14.output, sprintf("ribo_thres_%s", ribo.thres)),
                                                         use.sctransform = TRUE,
                                                         num.PCA = num.PCA,
                                                         num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                         num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                         cluster.resolution = cluster.resolution,
                                                         vars.to.regress = vars.to.regress,
                                                         k.weight = 50, 
                                                         PROJECT = PROJECT)  
    
  } else {
    print(sprintf(
      "Data exists at %s",
      file.path(path.to.14.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds")
    ))
    s.obj.integrated <- readRDS(file.path(path.to.14.output, 
                                          sprintf("ribo_thres_%s", ribo.thres), 
                                          "s8_output", 
                                          "integrated_BrainMet_dataset.output.s8.rds"))
  }
  
  #####----------------------------------------------------------------------#####
  ##### Update 13.01.2025
  #####----------------------------------------------------------------------#####
  path.to.main.src <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq/integrate_BrainMet_datasets"
  path.to.module.score <- file.path(path.to.main.src, "module_scores_20250113.xlsx")
  module.genedf <- readxl::read_excel(path.to.module.score)
  module.gene.list <- list()
  single.module.genes <- c()
  for (g in colnames(module.genedf)){
    tmp <- module.genedf[[g]] %>% unique()
    g <- str_replace_all(g, " ", "_")
    module.gene.list[[g]] <- tmp[is.na(tmp) == FALSE]
    single.module.genes <- c(single.module.genes, module.gene.list[[g]])
  }

  for (input.col in names(module.gene.list)){
    module.gene.list[[input.col]] <- unlist(lapply(module.gene.list[[input.col]], function(x){
      return(str_trim(x, "right"))
    }))
  }

  for (input.list in names(module.gene.list)){
    DefaultAssay(s.obj.integrated) <- "SCT"
    s.obj.integrated <- AddModuleScore(object = s.obj.integrated, features = list(module.gene.list[[input.list]]), name = sprintf("%s_", input.list), ctrl = 50)
  }

  fake.module.gene.list <- to_vec(
    for (item in names(module.gene.list)){
      sprintf("%s_1", item)
    }
  )

  DefaultAssay(s.obj.integrated) <- "SCT"
  Idents(s.obj.integrated) <- "harmony.cluster.0.5"

  dir.create(file.path(path.to.14.output, sprintf("ribo_thres_%s", ribo.thres), "module_scores"), showWarnings = FALSE, recursive = TRUE)

  for (input.module in fake.module.gene.list){

    feature.plot.module.scores <- FeaturePlot(object = s.obj.integrated,
                                              reduction = "harmony_UMAP",
                                              label = TRUE,
                                              features = c(input.module),
                                              pt.size = 0.5) &
      scale_color_gradient(low = "lightgray", high = "#FF0000", na.value = "lightgray")
    violin.plot.module.scores <- VlnPlot(object = s.obj.integrated,
                                         features = c(input.module),
                                         pt.size = 0,
                                         group.by = "harmony.cluster.0.5")
    ggsave(plot = feature.plot.module.scores + violin.plot.module.scores,
           filename = sprintf("UMAP_and_ViolinPlot_%s.svg", str_replace_all(input.module, "/", "_")),
           path = file.path(path.to.14.output, sprintf("ribo_thres_%s", ribo.thres), "module_scores"),
           width = 14,
           height = 10,
           dpi = 300,
           device = "svg")
  }

  for (input.module in fake.module.gene.list){
    dir.create(file.path(path.to.14.output, sprintf("ribo_thres_%s", ribo.thres), "single_genes", str_replace_all(input.module, "/", "_")),
               showWarnings = FALSE, recursive = TRUE)
    plot.genes <- intersect(row.names(s.obj.integrated), module.gene.list[[str_replace(input.module, "_1", "")]])
    for (input.gene in plot.genes){
      feature.plot.module.scores <- FeaturePlot(object = s.obj.integrated,
                                                reduction = "harmony_UMAP",
                                                label = TRUE,
                                                features = c(input.gene),
                                                pt.size = 0.5) &
        scale_color_gradient(low = "lightgray", high = "#FF0000", na.value = "lightgray")
      violin.plot.module.scores <- VlnPlot(object = s.obj.integrated,
                                           features = c(input.gene),
                                           pt.size = 0, group.by = "harmony.cluster.0.5")
      ggsave(plot = feature.plot.module.scores + violin.plot.module.scores,
             filename = sprintf("UMAP_and_ViolinPlot_%s.svg", str_replace_all(input.gene, "/", "_")),
             path = file.path(path.to.14.output, sprintf("ribo_thres_%s", ribo.thres), "single_genes", str_replace_all(input.module, "/", "_")),
             width = 14,
             height = 10,
             dpi = 300,
             device = "svg")
    }
  }
}

# removed.cells <- setdiff(colnames(s.obj.integrated.12), colnames(s.obj.integrated))
# meta.data <- s.obj.integrated.12@meta.data %>% rownames_to_column("barcode") %>% subset(barcode %in% removed.cells)
# 
# count.remove.cells <- table(meta.data$name) %>% sort()
# count.total.cells <- table(s.obj.integrated.12@meta.data$name)
# count.total.cells <- count.total.cells[names(count.remove.cells)]
# 
# VlnPlot(object = s.obj.integrated.12, features = c("percent.ribo"), group.by = "harmony.cluster.0.5")
