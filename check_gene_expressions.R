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

# view(dataset.metadata)
# row_i <- 25
row_i <- 31

info.col <- list()
for (c in colnames(dataset.metadata)){
  info.col[[c]] <- dataset.metadata[row_i, ][[c]]
}

info.col <- info.col[is.na(info.col) == FALSE]

save.name <- paste0(info.col[setdiff(names(info.col), c("path", "unique.name"))], collapse = "_" )
path.to.dge.output <- file.path(path.to.main.output2, "DGE", save.name)

path.to.s.obj <- info.col$path

print("reading in seurat object ...")
s.obj <- readRDS(str_replace(path.to.s.obj, ".rds", ".addLabel.rds"))

path.to.save.output <- file.path(path.to.main.output, "check_gene_expression_20250826", save.name)
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.save.output, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.save.output, "check_DE_genes"), showWarnings = FALSE, recursive = TRUE)

cluster.resolution <- 0.4
# sanity check: UMAP
p.umap <- DimPlot(object = s.obj, reduction = "harmony_UMAP", 
                  group.by = sprintf("harmony.cluster.%s", cluster.resolution),
                  label = TRUE, label.box = TRUE)
ggsave(plot = p.umap, filename = "UMAP.pdf", path = file.path(path.to.save.output, "figures"), device = "pdf", width = 14, height = 10, dpi = 300)

check.genes <- c("CD3E", "CD3D", "CD4", "CD8A")

Idents(s.obj) <- sprintf("harmony.cluster.%s", cluster.resolution)

for (input.name in unique(s.obj$label)){
  p.feature.plot <- FeaturePlot(object = subset(s.obj, label == input.name), 
                                reduction = "harmony_UMAP", 
                                features = check.genes, 
                                label = TRUE)
  ggsave(plot = p.feature.plot, 
         filename = sprintf("UMAP_genes_%s.pdf", input.name), 
         path = file.path(path.to.save.output, "figures"), 
         device = "pdf", 
         width = 14, 
         height = 10, 
         dpi = 300)
  
  p.violin.plot <- VlnPlot(object = subset(s.obj, label == input.name), features = check.genes, pt.size = 0)
  ggsave(plot = p.violin.plot, filename = sprintf("violin_plot_%s.pdf", input.name), 
         path = file.path(path.to.save.output, "figures"), 
         device = "pdf", 
         width = 14, 
         height = 10, 
         dpi = 300)
}



##### check why the genes AC... appear too many times in DGE

# for (item in row.names(s.obj)){
#   if (grepl("AC", item) == TRUE){
#     print(item)
#   }
# }


col <- "label"

all.conditions <- c("BrainMet",
                    "Control Epilepsy",
                    "Control Glioma")
all.comparisons <- data.frame(
  condition1 = c("BrainMet", "BrainMet", "Control Glioma"),
  condition2 = c("Control Epilepsy", "Control Glioma", "Control Epilepsy")
)


for (j in seq(1, nrow(all.comparisons))){
  sig.genes <- list()
  
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

  all.clusters <- Sys.glob(file.path(outputdir, "cluster_*")) %>% basename()
  
  all.sig.genes <- c()
  
  for (cluster.id in all.clusters){
    if (file.exists(file.path(outputdir, cluster.id, "sig_DE_gene_list.csv")) == TRUE){
      sig.genes[[cluster.id]] <- read.csv(file.path(outputdir, cluster.id, "sig_DE_gene_list.csv"))$gene     
      all.sig.genes <- c(all.sig.genes, sig.genes[[cluster.id]]) 
    }
  }
  countdf <- table(all.sig.genes) %>% data.frame() %>% arrange(desc(Freq))
  colnames(countdf) <- c("Gene", "count")
  writexl::write_xlsx(countdf, file.path(path.to.save.output, 
                                         sprintf("count_genes_in_DE_results_%s_vs_%s.xlsx", condition1, condition2)))
  
  # genes.in.all.de <- Reduce(intersect, sig.genes)
  
  p.compare <- VlnPlot(object = s.obj, features = head(countdf, 9)$Gene, group.by = "label", pt.size = 0, ncol = 3)
  ggsave(plot = p.compare, filename = sprintf("violin_plot_top9_frequent_genes_in_DE_%s_vs_%s.pdf", condition1, condition2),
         path = path.to.save.output, device = "pdf", width = 14, height = 10, dpi = 300)
}
  