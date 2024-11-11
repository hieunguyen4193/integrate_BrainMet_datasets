gc()
rm(list = ls())

options(future.globals.maxSize = 10000 * 1024^2)

if ("argparse" %in% installed.packages() == FALSE){
  install.packages("argparse")
}
library(argparse)
integrated.version <- "v0.1"
sample.id <- sprintf("integrated_%s", integrated.version)

my_random_seed <- 42
set.seed(my_random_seed)

# __________VDJ DATA ANYLYSIS PIPELINE__________
PROJECT <- "BrainMet_SeuratV5"
version.name <- sprintf("SeuratV5_%s", sample.id)

PROJECT.with.version <- sprintf("%s_%s", PROJECT, version.name)

outdir <- "/media/hieunguyen/HNSD_mini/data/outdir"

path.to.main.output <- file.path(outdir, PROJECT)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
path2src <- file.path(path.to.pipeline.src, "processes_src")

source(file.path(path2src, "s2_ambient_RNA_correction.R"))
source(file.path(path2src, "s3_first_filter.R"))
source(file.path(path2src, "s8_integration_and_clustering_SeuratV5.R"))

source(file.path(path2src, "import_libraries.R"))
source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))

analysis.round <- "1st"
cluster.resolution <- 0.5
cluster.resolution.small <- 0.02

path.to.main.output <- file.path(outdir, PROJECT, sprintf("integrated_%s", integrated.version))
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, sprintf("03_output_%s", cluster.resolution))
dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)

#####----------------------------------------------------------------------#####
##### S2: AmbientRNA estimate
#####----------------------------------------------------------------------#####
num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25
regressOut_mode <- NULL
features_to_regressOut <- NULL
use.sctransform <- FALSE
vars.to.regress <- c("percent.mt")
save.RDS.s8 <- TRUE
path.to.output <- path.to.03.output
all.reductions <- c("integrated.rpca")

if (file.exists(file.path(path.to.output, "s8_output", paste0(PROJECT, ".output.s8.rds"))) == FALSE){
  s.obj <- readRDS(file.path(path.to.02.output, "integrated.rds"))
  s.obj <- JoinLayers(s.obj)
  
  #####----------------------------------------------------------------------#####
  ##### S1: QC
  #####----------------------------------------------------------------------#####
  s.obj[["percent.mt"]] <- PercentageFeatureSet(s.obj, 
                                                pattern = "^mt-|^MT-")
  s.obj[["percent.ribo"]] <- PercentageFeatureSet(s.obj, 
                                                  pattern = "^Rpl|^Rps|^RPL|^RPS")
  
  s.obj <- s2_ambient_RNA_correction(s.obj = s.obj, 
                                     chosen.method = "decontX",
                                     path.to.output = path.to.03.output,
                                     save.RDS.s2 = TRUE,
                                     PROJECT = PROJECT)  
  
  #####----------------------------------------------------------------------#####
  ##### S3: filter
  #####----------------------------------------------------------------------#####
  filter.thresholds <- list( nFeatureRNAfloor = NULL,
                             nFeatureRNAceiling = NULL,
                             nCountRNAfloor = NULL,
                             nCountRNAceiling = NULL,
                             pct_mitofloor = NULL,
                             pct_mitoceiling = 10,
                             pct_ribofloor = NULL,
                             pct_riboceiling = NULL,
                             ambientRNA_thres = 0.5)
  
  s.obj <- s3.filter(s.obj = s.obj,
                     PROJECT = PROJECT,
                     path.to.output = path.to.03.output,
                     save.RDS.s3 = TRUE,
                     nFeatureRNAfloor = filter.thresholds$nFeatureRNAfloor, 
                     nFeatureRNAceiling = filter.thresholds$nFeatureRNAceiling,
                     nCountRNAfloor = filter.thresholds$nCountRNAfloor, 
                     nCountRNAceiling = filter.thresholds$nCountRNAceiling,
                     pct_mitofloor = filter.thresholds$pct_mitofloor, 
                     pct_mitoceiling = filter.thresholds$pct_mitoceiling,
                     pct_ribofloor = filter.thresholds$pct_ribofloor, 
                     pct_riboceiling = filter.thresholds$pct_riboceiling,
                     ambientRNA_thres = filter.thresholds$ambientRNA_thres,
                     log10GenesPerUMI_thres = filter.thresholds$log10GenesPerUMI)
  s.obj[["RNA"]] <- split(s.obj[["RNA"]], f = s.obj$name)
  DefaultAssay(s.obj) <- "RNA"
  
  saveRDS(s.obj, file.path(path.to.03.output, "s3_output", sprintf("%s.splitted.output.s3.rds", PROJECT)))
  
  if (use.sctransform == TRUE){
    print("Run SCTransform ...")
    s.obj <- SCTransform(s.obj, vars.to.regress = vars.to.regress, verbose = FALSE, conserve.memory = TRUE)
    normalization.method <- "SCT" 
  } else {
    print("Running normalization and scaling data")
    s.obj <- NormalizeData(s.obj, normalization.method = "LogNormalize") # ---> use Log Normalized
    s.obj <- FindVariableFeatures(s.obj, selection.method = "vst", nfeatures = 3000)
    s.obj <- ScaleData(s.obj, features = VariableFeatures(s.obj))
    normalization.method <- "LogNormalize" 
  }
  print("Running PCA and UMAP")
  s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = TRUE, reduction.name = "RNA_PCA")  
  s.obj <- RunUMAP(s.obj, dims = 1:num.PC.used.in.UMAP, reduction = "RNA_PCA", reduction.name = "umap.unintegrated")  
  
  saveRDS(s.obj, file.path(path.to.03.output, "s3_output", sprintf("%s.scaled.splitted.output.s3.rds", PROJECT)))
  
  if (use.sctransform == TRUE){
    normalization.method <- "SCT"
  } else {
    normalization.method <- "LogNormalize"
  }
  
  print("Integrate layers ...")
  
  
  print("start integration")
  
  s.obj <- IntegrateLayers(
    object = s.obj, 
    method = RPCAIntegration,
    orig.reduction = "RNA_PCA", 
    new.reduction = "integrated.rpca", 
    verbose = FALSE, 
    normalization.method = normalization.method
  )
  saveRDS(s.obj, file.path(path.to.03.output, "s3_output", sprintf("%s.pre_integrated.output.s3.rds", PROJECT)))
  
  for (selected.reduction in all.reductions){
    new.reduction.name <- str_replace(selected.reduction, "integrated.", "")
    print("RUnning Find neighbors")
    s.obj <- FindNeighbors(s.obj, dims = 1:num.PC.used.in.Clustering, reduction = selected.reduction)
    print("Running Find clusters")
    s.obj <- FindClusters(s.obj, resolution = cluster.resolution, cluster.name = sprintf("%s.cluster.%s", new.reduction.name, cluster.resolution))
    print("Running UMAP")
    s.obj <- RunUMAP(
      object = s.obj,
      reduction = selected.reduction,
      dims = 1:num.PC.used.in.UMAP, 
      reduction.name = sprintf("%s_UMAP", new.reduction.name))
  }
  
  print("Start saving ...")
  if (save.RDS.s8 == TRUE){
    dir.create(file.path(path.to.output, "s8_output"), showWarnings = FALSE, recursive = TRUE)
    saveRDS(object = s.obj, 
            file = file.path(path.to.output, "s8_output", 
                             paste0(PROJECT, ".output.s8.rds")))
  }
} else {
  print("reading in data at s3")
  s.obj <- readRDS(file.path(path.to.output, "s8_output", paste0(PROJECT, ".output.s8.rds")))
}

print(sprintf("Start running clustering with lower resolution at %s", cluster.resolution.small))
for (selected.reduction in all.reductions){
  s.obj <- FindNeighbors(s.obj, dims = 1:num.PC.used.in.Clustering, reduction = selected.reduction)
  s.obj <- FindClusters(s.obj, resolution = cluster.resolution.small, cluster.name = sprintf("%s.cluster.%s", new.reduction.name, cluster.resolution.small))
}

saveRDS(object = s.obj, file.path(path.to.output, "s8_output", paste0(PROJECT, "small_clusters.output.s8.rds")))

#####----------------------------------------------------------------------#####
##### INTEGRATION
#####----------------------------------------------------------------------#####

Idents(s.obj) <- "seurat_clusters"
if (file.exists(file.path(path.to.03.output, "DE_cluster_marker_genes.rds")) == FALSE){
  print("Start preparation for find all markers ...")
  DefaultAssay(s.obj) <- "RNA"
  s.obj <- JoinLayers(s.obj)
  print("start running FindAllMarkers ...")
  cluster.markers <- FindAllMarkers(object = subset(s.obj, seurat_clusters %in% seq(1,10)), assay = "RNA", test.use = "wilcox", slot = "data", min.pct = 0.5)
  cluster.markers <- subset(cluster.markers, cluster.markers$p_val_adj < 0.05 & cluster.markers$avg_log2FC > 0)
  saveRDS(cluster.markers, file.path(path.to.03.output, "DE_cluster_marker_genes.rds"))
} else {
  cluster.markers <- readRDS(file.path(path.to.03.output, "DE_cluster_marker_genes.rds"))
}