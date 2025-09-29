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

all.conditions <- c("BrainMet",
                    "Control Epilepsy",
                    "Control Glioma" )

if ("glmGamPoi" %in% installed.packages() == FALSE){
  BiocManager::install('glmGamPoi', update = FALSE)
}

for (input.condition in all.conditions){
  path.to.12.output <- file.path(path.to.main.output, "12_output_separated_conditions", input.condition)
  dir.create(path.to.12.output, showWarnings = FALSE, recursive = TRUE)
  
  if (file.exists(file.path(path.to.12.output, "s8_output", sprintf("%s.output.s8.rds", input.condition))) == FALSE){
    # read data from the v0.2 dataset and the GSE193745 dataset
    s.obj.GSE193745 <- readRDS(file.path(path.to.11.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
    s.obj.prev <- readRDS(file.path(path.to.09.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
    
    selected.cells <- c(colnames(s.obj.GSE193745), colnames(s.obj.prev))
    
    ##### get data from v0.3 integration version
    s.obj <- readRDS(file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", "v0.3"), "02_output", "integrated.rds"))
    
    s.obj <- subset(s.obj, cells = selected.cells)
    print(sprintf("Number of cells left after filtering: %s", length(colnames(s.obj))))
    
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
    
    # add conditions (Brainmet, control) to the seurat object
    sample.metadata <- read.csv(file.path(path.to.main.src, "samples_to_integrated_20240725.csv"))
    meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
    meta.data <- merge(meta.data, sample.metadata, by.x = "name", by.y = "Sample")
    meta.data <- meta.data %>% column_to_rownames("barcode")
    meta.data <- meta.data[row.names(s.obj@meta.data), ]
    s.obj <- AddMetaData(object = s.obj, col.name = "label", metadata = meta.data$Group)
    
    # subseting the seurat object based on conditions
    print(sprintf("subsetting the seurat object to condition %s", input.condition))
    s.obj <- subset(s.obj, label == input.condition)
    
    DefaultAssay(s.obj) <- "RNA"
    print("Running join layers ...")
    s.obj <- JoinLayers(s.obj)
    
    s.obj.integrated <- s8.integration.and.clustering_V5(s.obj = s.obj, 
                                                         save.RDS.s8 = TRUE,
                                                         path.to.output = path.to.12.output,
                                                         use.sctransform = TRUE,
                                                         num.PCA = num.PCA,
                                                         num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                         num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                         cluster.resolution = cluster.resolution,
                                                         vars.to.regress = vars.to.regress,
                                                         PROJECT = input.condition)  
  } else {
    print("reading in saved seurat object")
    s.obj.integrated <- readRDS(file.path(path.to.12.output, "s8_output", sprintf("%s.output.s8.rds", input.condition)))
    print("finished saving data ...")
  }
  
  s.obj <- s.obj.integrated
  #####--------------------------------------------------------------------#####
  ##### find cluster markers and generate more cluster rsolution data
  #####--------------------------------------------------------------------#####
  print("start running cluster markers and more cluster resolution clustering ...")
  
  selected.reduction <- "harmony"
  num.PCA <- 25
  num.PC.used.in.UMAP <- 25
  num.PC.used.in.Clustering <- 25
  regressOut_mode <- NULL
  features_to_regressOut <- NULL
  use.sctransform <- TRUE
  vars.to.regress <- c("percent.mt")
  
  path.to.save.output <- path.to.12.output
  if (file.exists(file.path(path.to.save.output, "integrated_BrainMet_dataset.moreClusterRes.output.s8.rds")) == FALSE){
    for (cluster.resolution in c(0.1, 0.2, 0.25, 0.3, 0.4)){
      print("-----------------------------------------------------------")
      print(sprintf("Working on clustering with resolution %s", cluster.resolution))
      print("-----------------------------------------------------------")
      s.obj <- FindNeighbors(s.obj, dims = 1:num.PC.used.in.Clustering, reduction = selected.reduction)
      s.obj <- FindClusters(s.obj, resolution = cluster.resolution, cluster.name = sprintf("%s.cluster.%s", selected.reduction, cluster.resolution))  
    }
    saveRDS(s.obj, file.path(path.to.save.output, "integrated_BrainMet_dataset.moreClusterRes.output.s8.rds")) 
  } else {
    print("reading in saved data ...")
    s.obj <- readRDS(file.path(path.to.save.output, "integrated_BrainMet_dataset.moreClusterRes.output.s8.rds"))
  }
  
  for (cluster.resolution in c(0.1, 0.2, 0.25, 0.3, 0.4)){
    dir.create(file.path(path.to.save.output, sprintf("cluster_resolution_%s", cluster.resolution)), showWarnings = FALSE, recursive = TRUE)
    ##### save all UMAP cor. all cluster resolutions
    p <- DimPlot(object = s.obj, reduction = "harmony_UMAP", label = TRUE, label.box = TRUE, 
                 group.by = sprintf("harmony.cluster.%s", cluster.resolution))
    ggsave(plot = p, filename = sprintf("12a_output_integrated_BrainMet.harmony.cluster.%s.pdf", cluster.resolution),
           path = file.path(path.to.save.output, sprintf("cluster_resolution_%s", cluster.resolution)), device = "pdf", width = 14, height = 10, dpi =  300)    
    
    Idents(s.obj) <- sprintf("harmony.cluster.%s", cluster.resolution)
    min.pct <- 0.1
    s.obj <- PrepSCTFindMarkers(s.obj)
    print("start running FindAllMarkers ...")
    cluster.markers <- FindAllMarkers(object = s.obj, 
                                      assay = "SCT", 
                                      test.use = "wilcox", 
                                      slot = "data", 
                                      min.pct = min.pct, 
                                      recorrect_umi = TRUE)
    cluster.markers <- subset(cluster.markers, 
                              cluster.markers$p_val_adj < 0.05 & 
                                cluster.markers$avg_log2FC > 0)
    for (cluster.id in unique(cluster.markers$cluster)){
      writexl::write_xlsx(
        subset(cluster.markers, cluster.markers$cluster == cluster.id) %>% arrange(desc(avg_log2FC)),
        file.path(path.to.save.output, sprintf("cluster_resolution_%s", cluster.resolution), sprintf("cluster_%s.xlsx", cluster.id))
      )
    }
  }
}





