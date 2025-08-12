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

path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version), "scAnnotationPipeline")
path.to.01.output <- file.path(path.to.main.output, "outdir")
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

path.to.s.obj <- "/media/hieunguyen/GSHD_HN01/outdir/BrainMet_SeuratV5/integrate_BrainMet_datasets/integrated_v0.2/12a_output/integrated_BrainMet_dataset.moreClusterRes.output.s8.rds"

s.obj <- readRDS(path.to.s.obj)

all.cluster.names <- to_vec(
  for (item in colnames(s.obj@meta.data)){
    if (grepl("harmony.cluster", item) == TRUE){
      item
    }
  }
)

#####----------------------------------------------------------------------#####
##### cell type annotation by find all clusters' markers
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.01.output, "check_FindAllMarkers.csv")) == FALSE){
  min.pct <- 0.5
  for (selected.cluster in all.cluster.names){
    DefaultAssay(s.obj) <- "SCT"
    s.obj <- PrepSCTFindMarkers(s.obj)
    print("start running FindAllMarkers ...")
    Idents(s.obj) <- selected.cluster
    cluster.markers <- FindAllMarkers(object = s.obj, assay = "SCT", test.use = "wilcox", slot = "data", min.pct = min.pct, recorrect_umi = TRUE)
    cluster.markers.raw <- cluster.markers
    cluster.markers <- subset(cluster.markers, cluster.markers$p_val_adj < 0.05 & cluster.markers$avg_log2FC > 0)
    saveRDS(cluster.markers, file.path(path.to.01.output, sprintf("cluster_markers_cluster_res_%s.rds", selected.cluster)))  
    saveRDS(cluster.markers.raw, file.path(path.to.01.output, sprintf("cluster_markers_cluster_res_%s.raw.rds", selected.cluster)))  
  }
  write.csv(data.frame(status = c("check_finished_FindAllMarkers.csv")), file.path(path.to.01.output, "check_FindAllMarkers.csv"))
} 

#####----------------------------------------------------------------------#####
##### annotation with SingleR and Celldex
#####----------------------------------------------------------------------#####
if ("SingleR" %in% installed.packages() == FALSE){
  BiocManager::install("SingleR")
} 
if ("celldex" %in% installed.packages() == FALSE){
  BiocManager::install("celldex")
  # install the newest celldex database, v1.18, not 1.12
  # BiocManager::install(c("alabaster.base", "alabaster.matrix", "alabaster.se"))
  # install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/gypsum_1.4.0.tar.gz", type = "source", repos = NULL)
  # install.packages("https://bioconductor.org/packages/release/data/experiment/src/contrib/celldex_1.18.0.tar.gz",type = "source", repos = NULL)
} 

library(SingleR)
library(celldex)

db <- celldex::DatabaseImmuneCellExpressionData()
singleR.preddf <- SingleR(test = as.SingleCellExperiment(s.obj), ref = db, assay.type.test=1,
                     labels = db$label.main)
singleR.preddf <- data.frame(singleR.preddf)
write.csv(singleR.preddf, file.path(path.to.01.output, "SingleR_Celldex_ImmuneCell_prediction.csv"))

