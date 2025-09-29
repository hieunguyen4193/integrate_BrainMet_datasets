gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### data preparation
#####----------------------------------------------------------------------#####

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))
library(bluster)
path.to.main.src <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq/integrate_BrainMet_datasets"

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

evaluation.version <- "20250928"
# recently added, keep many versions. 
path.to.save.output <- file.path(path.to.main.output, "integration_evaluation", evaluation.version)
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

sep.annotations <- readxl::read_xlsx(file.path(path.to.main.src, "12SEP_cell_annotations_0.4.xlsx"))

all.s.obj <- list()

all.conditions <- c("BrainMet",
                    "Control Epilepsy",
                    "Control Glioma" )

##### reading in saved data
# integrated data 
s.obj.integrated <- readRDS(file.path(path.to.main.output, "12a_output", "integrated_BrainMet_dataset.moreClusterRes.output.s8.rds"))

# separated integrated data 
for (input.condition in all.conditions){
  print(sprintf("reading in input condition %s", input.condition))
  path.to.12.output <- file.path(path.to.main.output, "12_output_separated_conditions", input.condition)
  tmp.s.obj <- readRDS(file.path(path.to.12.output, "integrated_BrainMet_dataset.moreClusterRes.output.s8.rds"))
  tmp.annotdf <- subset(sep.annotations, sep.annotations$label  == input.condition)
  
  tmp.metadata <- tmp.s.obj@meta.data %>% rownames_to_column("barcode") %>%
    rowwise() %>%
    mutate(celltype_and_origin = sprintf("%s_%s", input.condition, subset(tmp.annotdf, 
                                                                            tmp.annotdf$cluster == harmony.cluster.0.4 & 
                                                                            tmp.annotdf$label == input.condition)$celltype)) %>%
    mutate(celltype = subset(tmp.annotdf, tmp.annotdf$cluster == harmony.cluster.0.4 )$celltype) %>%
    column_to_rownames("barcode")
  tmp.metadata <- tmp.metadata[row.names(tmp.s.obj@meta.data), ]
  
  tmp.s.obj <- AddMetaData(object = tmp.s.obj, col.name = "celltype", metadata = tmp.metadata$celltype)
  tmp.s.obj <- AddMetaData(object = tmp.s.obj, col.name = "celltype_and_origin", metadata = tmp.metadata$celltype_and_origin)
  all.s.obj[[input.condition]]  <- tmp.s.obj
}

merge.metadata <- data.frame()
for (input.condition in all.conditions){
  tmp.metadata <- all.s.obj[[input.condition]]@meta.data %>% rownames_to_column("barcode") %>%
    subset(select = c(barcode, celltype))
  merge.metadata <- rbind(merge.metadata, tmp.metadata)
}

integrated.metadata <- s.obj.integrated@meta.data %>% rownames_to_column("barcode")

integrated.metadata <- merge(integrated.metadata, merge.metadata, by.x = "barcode", by.y = "barcode") %>% column_to_rownames("barcode")
integrated.metadata <- integrated.metadata[row.names(s.obj.integrated@meta.data), ]

s.obj.integrated <- AddMetaData(object = s.obj.integrated, metadata = integrated.metadata$celltype, col.name = "celltype")

metadata.barcode.celltype <- s.obj.integrated@meta.data %>% rownames_to_column("barcode")
write.csv(metadata.barcode.celltype, file.path(path.to.save.output, "metadata_barcode_celltype.csv"))

if ("ggthemes" %in% installed.packages() == FALSE){
  install.packages("ggthemes")
}
library(ggthemes)
input.colors <- tableau_color_pal(
  palette = "Tableau 20",
  type = c("regular", "ordered-sequential", "ordered-diverging"),
  direction = 1
)

group.by.input <- "celltype"

# plot UMAP with cell type from sep integration
p.umap.integrated <- DimPlot(object = s.obj.integrated, label = TRUE, label.box = TRUE, group.by = group.by.input, reduction = "harmony_UMAP", pt.size = 1, 
                             cols = input.colors(length(unique(s.obj.integrated@meta.data[[group.by.input]]))) ) +
  xlim(c(-10, 15)) + ylim(c(-8, 8))

meta.data <- s.obj.integrated@meta.data %>% rownames_to_column("barcode")

countdf <- table(meta.data$harmony.cluster.0.4, meta.data$celltype) %>% data.frame()
colnames(countdf) <- c("cluster", "label", "count")
sum.count <- list()
for (i in unique(countdf$cluster)){
  sum.count[[sprintf("cluster_%s", i)]] <- subset(countdf, countdf$cluster == i)$count %>% sum()
}
countdf <- countdf %>% rowwise() %>%
  mutate(pct = count/sum.count[[sprintf("cluster_%s", cluster)]])

p.count <- countdf %>% ggplot(aes(x = cluster, y = count, fill = label)) + geom_bar(stat = "identity")

# composition of cell types in each cluster. 
p.pct <- countdf %>% ggplot(aes(x = cluster, y = pct, fill = label)) + geom_bar(stat = "identity") + 
  scale_fill_manual( values = input.colors(length(unique(s.obj.integrated@meta.data[[group.by.input]])))) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

#####----------------------------------------------------------------------#####
##### Calculate metrics to evaluate the integration
#####----------------------------------------------------------------------#####
sample.metadata <- read.csv(file.path(path.to.project.src, "samples_to_integrated_20240725.csv"))
meta.data <- s.obj.integrated@meta.data %>% rownames_to_column("barcode")
meta.data <- merge(meta.data, sample.metadata, by.x = "name", by.y = "Sample")
meta.data <- meta.data %>% column_to_rownames("barcode")
meta.data <- meta.data[row.names(s.obj.integrated@meta.data), ]

s.obj.integrated <- AddMetaData(object = s.obj.integrated, col.name = "label", metadata = meta.data$Group)

# plot UMAP with same cell type from different conditions
dir.create(file.path(path.to.save.output, "UMAP"), showWarnings = FALSE, recursive = TRUE)

for (input.celltype in unique(s.obj.integrated$celltype)){
  p.umap.celltype <- DimPlot(object = subset(s.obj.integrated, celltype == input.celltype), 
                             label = TRUE, 
                             label.box = TRUE, 
                             group.by = "label", 
                             reduction = "harmony_UMAP", 
                             pt.size = 1) +
    xlim(c(-10, 15)) + ylim(c(-8, 8)) +
    ggtitle(input.celltype)  
  ggsave(plot = p.umap.celltype, filename = sprintf("%s.pdf", input.celltype), 
         path = file.path(path.to.save.output, "UMAP"), device = "pdf", width = 14, height = 10, dpi = 300)
}

# input.condition <- "BrainMet"
# input.condition <- "Control Epilepsy" 
# input.condition <- "Control Glioma"
for (input.condition in all.conditions){
  for (selected.cluster in c("harmony.cluster.0.4", "celltype")){
    tmp.s.obj.integrated <- subset(s.obj.integrated, label == input.condition)
    tmp.s.obj.old <- all.s.obj[[input.condition]]
    
    meta.data1 <- tmp.s.obj.integrated@meta.data %>% rownames_to_column("barcode")
    meta.data1 <-  meta.data1[, c("barcode", "harmony.cluster.0.4")]
    colnames(meta.data1) <- c("barcode", "new.cluster")
    
    meta.data2 <- tmp.s.obj.old@meta.data %>% rownames_to_column("barcode") 
    meta.data2 <-  meta.data2[, c("barcode", selected.cluster)]
    
    colnames(meta.data2) <- c("barcode", "old.cluster")
    
    check.metadata <- merge(meta.data1, meta.data2, by.x = "barcode", by.y = "barcode")
    
    compare.old.new.clusters <- table(check.metadata$old.cluster, check.metadata$new.cluster) %>% 
      data.frame() %>%
      pivot_wider(names_from = "Var1", values_from = "Freq")
    
    cluster_labels_b1_before <- check.metadata$old.cluster
    
    cluster_labels_b1_after <- check.metadata$new.cluster
    
    b1_ari <- pairwiseRand(
      cluster_labels_b1_before,
      cluster_labels_b1_after,
      mode = "index",
      adjusted = TRUE
    )
    
    print(sprintf("Biological heterogeneity conservation ARI in %s: %s",
                  input.condition, b1_ari))
    
    tab_b1 <- nestedClusters(ref = paste("Before", cluster_labels_b1_before),
                             alt = paste("After", cluster_labels_b1_after))
    
    pheatmap::pheatmap(tab_b1$proportions, cluster_row = FALSE,
                       cluster_col = FALSE, main = sprintf("%s, compare before-after clusters, SUM ROW = 1, ARI = %s", input.condition, round(b1_ari, 2)),
                       silent = TRUE, display_numbers = TRUE, 
                       filename = file.path(path.to.save.output, sprintf("heatmap_%s_before_%s_vs_after_clustering.pdf", selected.cluster, input.condition)))
    
    tab_b2 <- nestedClusters(alt = paste("Before", cluster_labels_b1_before),
                             ref = paste("After", cluster_labels_b1_after))
    
    pheatmap::pheatmap(tab_b2$proportions %>% t(), cluster_row = FALSE,
                       cluster_col = FALSE, main = sprintf("%s, compare before-after clusters, SUM COL = 1, ARI = %s", input.condition, round(b1_ari, 2)),
                       silent = TRUE, display_numbers = TRUE,
                       filename = file.path(path.to.save.output, sprintf("heatmap_%s_before_%s_vs_after_clustering.reverse.pdf", selected.cluster, input.condition)))
    
  }
}

# batch mixing 
batch_labels <- s.obj.integrated@meta.data$name
integrated_cluster_labels <- s.obj.integrated@meta.data$harmony.cluster.0.4
batch_ari <- 1 - pairwiseRand(batch_labels, integrated_cluster_labels,
                              mode = "index", adjusted = TRUE)
print(paste0("Batch mixing ARI: ", batch_ari))

batch_labels <- s.obj.integrated@meta.data$label
integrated_cluster_labels <- s.obj.integrated@meta.data$harmony.cluster.0.4
batch_ari <- 1 - pairwiseRand(batch_labels, integrated_cluster_labels,
                              mode = "index", adjusted = TRUE)
print(paste0("Batch mixing ARI: ", batch_ari))

# LISI: local inverse simpson index.
if ("lisi" %in% installed.packages() == FALSE){
  devtools::install_github("immunogenomics/lisi")
}

library(lisi)

umapdf <- s.obj.integrated@reductions$harmony_UMAP@cell.embeddings %>%
  data.frame() %>% rownames_to_column("barcode") 

umapdf <- merge(umapdf, s.obj.integrated@meta.data %>% rownames_to_column("barcode"))
lisi.metadata <- umapdf[, c("barcode", "celltype", "name")] %>%
  column_to_rownames("barcode")

lisidf <- compute_lisi(umapdf[, c("harmonyUMAP_1", "harmonyUMAP_2")], 
                       lisi.metadata, 
                       c("celltype", "name")) %>%
  rownames_to_column("barcode")
colnames(lisidf) <- c("barcode", "lisi_celltype", "lisi_name")

meta.data <- s.obj.integrated@meta.data %>% rownames_to_column("barcode")
meta.data <- merge(meta.data, lisidf, by.x = "barcode", by.y = "barcode") %>%
  column_to_rownames("barcode")

meta.data <- meta.data[row.names(s.obj.integrated@meta.data), ]

s.obj.integrated <- AddMetaData(object = s.obj.integrated, col.name = "lisi_celltype", metadata = meta.data$lisi_celltype)
s.obj.integrated <- AddMetaData(object = s.obj.integrated, col.name = "lisi_name", metadata = meta.data$lisi_name)

p <- FeaturePlot(object = s.obj.integrated, reduction = "harmony_UMAP", features = c("lisi_name"), label = TRUE)
ggsave(plot = p, filename = "lisi_name.pdf", path = path.to.save.output, device = "pdf", width = 14, height = 10)

p <- VlnPlot(object = s.obj.integrated, features = c("lisi_name"), pt.size = 0)
ggsave(plot = p, filename = "lisi_name.violinplot.pdf", path = path.to.save.output, device = "pdf", width = 14, height = 10)

p <- FeaturePlot(object = s.obj.integrated, reduction = "harmony_UMAP", features = c("lisi_celltype"), label = TRUE)
ggsave(plot = p, filename = "lisi_celltype.pdf", path = path.to.save.output, device = "pdf", width = 14, height = 10)

p <- VlnPlot(object = s.obj.integrated, features = c("lisi_celltype"), pt.size = 0)
ggsave(plot = p, filename = "lisi_celltype.violinplot.pdf", path = path.to.save.output, device = "pdf", width = 14, height = 10)


Idents(s.obj.integrated) <- "celltype"
min.pct <- 0.1 # raised to 0.1 on 28.09.2025

##### find cluster markers again
if (file.exists(file.path(path.to.save.output, "celltype_markers.rds")) == FALSE){
  cluster.markers <-   cluster.markers <- FindAllMarkers(object = s.obj.integrated, 
                                                         assay = "SCT", 
                                                         test.use = "wilcox", 
                                                         slot = "data", 
                                                         min.pct = min.pct, 
                                                         recorrect_umi = FALSE)
  saveRDS(cluster.markers, file.path(path.to.save.output, "celltype_markers.rds"))  
} else {
  print("cluster markers (celltype) exists, reading in again ...")
  cluster.markers <- readRDS(file.path(path.to.save.output, "celltype_markers.rds"))
}


dir.create(file.path(path.to.save.output, "cluster_markers"), showWarnings = FALSE, recursive = TRUE)

##### check CD3 genes
p <- FeaturePlot(s.obj.integrated, features = c("CD3D", "CD3E", "CD8A", "CD4"), label = TRUE, ncol = 2, reduction = "harmony_UMAP")
ggsave(plot = p, filename = "check_CD3_genes.pdf", path = path.to.save.output, device = "pdf", width = 14, height = 10)

##### cluster marker genes
for (cluster.id in unique(cluster.markers$cluster)){
  tmpdf <- subset(cluster.markers, cluster.markers$cluster == cluster.id) %>%
    subset(p_val_adj <= 0.05 & avg_log2FC >= 0) %>%
    arrange(desc(avg_log2FC))
  writexl::write_xlsx(tmpdf, file.path(path.to.save.output, "cluster_markers", sprintf("cluster_%s.xlsx", cluster.id)))
}

cell.countdf <- table(s.obj.integrated$label, s.obj.integrated$celltype) %>% data.frame() %>%
  pivot_wider(names_from = "Var1", values_from = "Freq")
writexl::write_xlsx(cell.countdf, file.path(path.to.save.output, "cell_count_per_conditions.xlsx"))

p <- DimPlot(object = s.obj.integrated, reduction = "harmony_UMAP", label = TRUE, label.box = TRUE, group.by = "celltype")
ggsave(plot = p, device = "pdf", width = 14, height = 10, dpi = 300, path = path.to.save.output, filename = "UMAP.celltype.pdf")

##### EOF