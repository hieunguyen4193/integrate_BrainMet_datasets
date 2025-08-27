gc()
rm(list = ls())

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

output.version <- "20240708"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq/integrate_BrainMet_datasets"
path.to.rmd <- file.path(path.to.project.src, "04_downstream_analysis_separate_conditions.Rmd")

all.integrated.metadata <- list(v0.1 = read.csv(file.path(path.to.project.src, "samples_to_integrated_20240513.csv")),
                                v0.2 = read.csv(file.path(path.to.project.src, "samples_to_integrated_20240620.csv")))
integrated.version <- "v0.2"
# outdir <- "/media/hieunguyen/HNSD01/outdir"
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
  path.to.17.output <- file.path(path.to.main.output2, "17_output", save.name)
  
  dir.create(path.to.17.output, showWarnings = FALSE, recursive = TRUE)
  if (file.exists(file.path(path.to.17.output, "cellchat_analysis.html")) == FALSE){
  # reading main seurat object. 
  s.obj <- readRDS(info.col$path)
  
  cluster.resolution <- info.col$cluster_resolution
  cluster.name <- sprintf("harmony.cluster.%s", cluster.resolution)
  reduction.name <- "harmony_UMAP"
  
  # plot the umap again for sanity check
  p.umap <- DimPlot(object = s.obj, reduction = reduction.name, group.by = cluster.name, label = TRUE, label.box = TRUE)
  ggsave(plot = p.umap, filename = sprintf("UMAP.pdf"), path = path.to.17.output, device = "pdf", width = 14, height = 10)
  
  # check metadata
  sample.metadata <- read.csv(file.path(path.to.project.src, "samples_to_integrated_20240725.csv"))
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
  meta.data <- merge(meta.data, sample.metadata, by.x = "name", by.y = "Sample")
  meta.data <- meta.data %>% column_to_rownames("barcode")
  meta.data <- meta.data[row.names(s.obj@meta.data), ]
  
  s.obj <- AddMetaData(object = s.obj, col.name = "label", metadata = meta.data$Group)
  
  p.umap.split.conditions <- DimPlot(object = s.obj, reduction = "harmony_UMAP", split.by = "label", ncol = 2, label = TRUE, label.box = TRUE)
  ggsave(plot = p.umap.split.conditions, filename = sprintf("UMAP_split_conditions.pdf"), path = path.to.17.output, device = "pdf", width = 14, height = 10)
  
  cell.countdf <- table(s.obj$label, s.obj@meta.data[[cluster.name]]) %>%
    data.frame() %>%
    pivot_wider(names_from = "Var1", values_from = "Freq")
  colnames(cell.countdf) <- c(cluster.name, "BrainMet", "Control Epilepsy", "Control Glioma")
  writexl::write_xlsx(cell.countdf, file.path(path.to.17.output, "cell_count_clusters_conditions.xlsx"))
  
  # Run cellchat analysis
  path.to.rmd <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq/integrate_BrainMet_datasets/cellchat_analysis_template.Rmd"
  
  
    rmarkdown::render(input = path.to.rmd, 
                      params = list(
                        reduction.name = reduction.name,
                        cluster.name = cluster.name,
                        path.to.s.obj = info.col$path,
                        outputdir = path.to.17.output
                      ),
                      output_file = "cellchat_analysis.html",
                      output_dir = path.to.17.output)    
  } else {
    print("data exists")
  }
}

# EOF