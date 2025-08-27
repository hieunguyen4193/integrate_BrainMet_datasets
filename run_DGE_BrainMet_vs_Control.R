gc()
rm(list = ls())

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
sctoolbox <- "/media/hieunguyen/HNSD01/src/sctoolbox" 

source(file.path(sctoolbox, "src/DGE/run_DGE_clusterwise.R"))
source(file.path(sctoolbox, "src/DGE/run_pathway_analysis.R"))
source(file.path(sctoolbox, "src/DGE/viz_pathway_analysis.R"))

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

output.version <- "20240708"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq/integrate_BrainMet_datasets"
path.to.rmd <- file.path(path.to.project.src, "04_downstream_analysis_separate_conditions.Rmd")

all.integrated.metadata <- list(v0.1 = read.csv(file.path(path.to.project.src, "samples_to_integrated_20240513.csv")),
                                v0.2 = read.csv(file.path(path.to.project.src, "samples_to_integrated_20240620.csv")))
integrated.version <- "v0.2"

outdir <- "/media/hieunguyen/GSHD_HN01/outdir"
outdir2 <- "/media/hieunguyen/HNSD01/outdir"
main.PROJECT <- "BrainMet_SeuratV5"

sample.metadata <- all.integrated.metadata[[integrated.version]]

main.PROJECT <- "BrainMet_SeuratV5"
code.version <- "integrate_BrainMet_datasets"
PROJECT <- "integrated_BrainMet_dataset"

path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))
path.to.main.output2 <- file.path(outdir2, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))

output.version <- '20240730'
path.to.save.html <- file.path(path.to.main.output, "html_outputs", output.version)
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.03.output <- file.path(path.to.main.output, "03_output")
path.to.09.output <- file.path(path.to.main.output, "09_output")
path.to.10.output <- file.path(path.to.main.output, "10_output")
path.to.11.output <- file.path(path.to.main.output, "11_output")
path.to.12.output <- file.path(path.to.main.output, "12_output")
path.to.13.output <- file.path(path.to.main.output, "13_output")
path.to.14.output <- file.path(path.to.main.output, "14_output")

dataset.metadata <- readxl::read_excel(file.path(path.to.project.src, "dataset_metadata.xlsx")) %>%
  rowwise() %>%
  mutate(unique.name = sprintf("%s_%s_%s_%s_%s",
                               output_index,
                               cluster_resolution,
                               ribo_filter,
                               dataset_name,
                               path))


for (row_i in seq(1, nrow(dataset.metadata))){
  info.col <- list()
  for (c in colnames(dataset.metadata)){
    info.col[[c]] <- dataset.metadata[row_i, ][[c]]
  }
  
  info.col <- info.col[is.na(info.col) == FALSE]
  
  save.name <- paste0(info.col[setdiff(names(info.col), c("path", "unique.name"))], collapse = "_" )
  path.to.dge.output <- file.path(path.to.main.output2, "DGE", save.name)
  
  dir.create(path.to.dge.output, showWarnings = FALSE, recursive = TRUE)
  
  path.to.s.obj <- info.col$path
  
  if (file.exists(str_replace(path.to.s.obj, ".rds", ".addLabel.rds")) == FALSE){
    s.obj <- readRDS(path.to.s.obj)
    sample.metadata <- read.csv(file.path(path.to.project.src, "samples_to_integrated_20240725.csv"))
    meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
    meta.data <- merge(meta.data, sample.metadata, by.x = "name", by.y = "Sample")
    meta.data <- meta.data %>% column_to_rownames("barcode")
    meta.data <- meta.data[row.names(s.obj@meta.data), ]
    s.obj <- AddMetaData(object = s.obj, col.name = "label", metadata = meta.data$Group)
    saveRDS(s.obj, str_replace(path.to.s.obj, ".rds", ".addLabel.rds"))    
  } else {
    print("file exists")
    print(str_replace(path.to.s.obj, ".rds", ".addLabel.rds"))
  }
  
  col <- "label"
  
  all.conditions <- c("BrainMet",
                      "Control Epilepsy",
                      "Control Glioma")
  all.comparisons <- data.frame(
    condition1 = c("BrainMet", "BrainMet", "Control Glioma"),
    condition2 = c("Control Epilepsy", "Control Glioma", "Control Epilepsy")
  )
  
  for (j in seq(1, nrow(all.comparisons))){
    condition1 <- all.comparisons[j, ]$condition1
    condition2 <- all.comparisons[j, ]$condition2
    
    cluster.resolution <- info.col$cluster_resolution
    cluster.name <- sprintf("harmony.cluster.%s", cluster.resolution)
    reduction.name <- "harmony_UMAP"
    cluster.id <- "all"
    remove.genes <- NULL
    min.pct <- 0
    assay <- "SCT"
    test.use <- "wilcox"
    
    outputdir <- file.path(
      path.to.dge.output, 
      sprintf("%s_vs_%s", condition1, condition2)
    )
    dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
    
    if (file.exists(file.path(outputdir, "SessionInfo.txt")) == FALSE){
      print(sprintf("Working on %s", save.name))
      print(sprintf("Working on sample %s vs sample %s", condition1, condition2))
      
      run_DGE_clusterwise(path.to.s.obj = str_replace(path.to.s.obj, ".rds", ".addLabel.rds"),
                          outputdir = outputdir,
                          condition1 = condition1,
                          condition2 = condition2,
                          col = col,
                          assay = "SCT",
                          reduction.name = reduction.name,
                          test.use = "wilcox",
                          cluster.id = cluster.id, 
                          cluster.name = cluster.name,
                          p_value_cutoff = 0.05,
                          min.pct = min.pct,
                          remove.genes = remove.genes)      
    } else {
      print("Data exists!")
    }
  }  
}




