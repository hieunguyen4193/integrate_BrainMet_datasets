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

sep.annotations <- readxl::read_xlsx(file.path(path.to.main.src, "12SEP_cell_annotations_0.4.xlsx"))

all.s.obj <- list()

all.conditions <- c("BrainMet",
                    "Control Epilepsy",
                    "Control Glioma" )

##### reading in saved data
# integrated data 
s.obj.integrated <- readRDS(file.path(path.to.main.output, "12a_output", "integrated_BrainMet_dataset.moreClusterRes.output.s8.rds"))

path.to.save.output <- file.path(path.to.main.output, "12_output_separated_conditions", "selected_cells_DGE")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

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

##### add label
path.to.project.src <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq/integrate_BrainMet_datasets"
sample.metadata <- read.csv(file.path(path.to.project.src, "samples_to_integrated_20240725.csv"))
meta.data <- s.obj.integrated@meta.data %>% rownames_to_column("barcode")
meta.data <- merge(meta.data, sample.metadata, by.x = "name", by.y = "Sample")
meta.data <- meta.data %>% column_to_rownames("barcode")
meta.data <- meta.data[row.names(s.obj.integrated@meta.data), ]

s.obj.integrated <- AddMetaData(object = s.obj.integrated, col.name = "label", metadata = meta.data$Group)

# finished loading data. 

##### MAIN ANALYSIS
selected.label <- "Control Glioma"
selected.celltype <- "Microglia"

dir.create(file.path(path.to.save.output, selected.label, selected.celltype), showWarnings = FALSE, recursive = TRUE)

s.obj <- subset(s.obj.integrated, celltype == selected.celltype) 
s.obj <- subset(s.obj, cells = colnames(all.s.obj[[selected.label]]))

p <- DimPlot(object = s.obj, reduction = "harmony_UMAP", label = TRUE, label.box = TRUE, group.by = "harmony.cluster.0.4")
ggsave(plot = p, filename = sprintf("%s_%s.pdf", selected.celltype, selected.label), 
       path = file.path(path.to.save.output, selected.label, selected.celltype), 
       device = "pdf", 
       width = 14, 
       height = 10, 
       dpi = 300)

if (selected.celltype == "TAMs"){
  cluster1 <- c(1)
  cluster2 <- c(2,5,7,9,14)  
} else if (selected.celltype == "Microglia"){
  cluster1 <- c(1, 4, 16)
  cluster2 <- c(2, 5, 7, 9, 14)
}

s.obj.subset <- subset(s.obj, harmony.cluster.0.4 %in% c(cluster1, cluster2))

tmp.metadata <- s.obj.subset@meta.data %>% rownames_to_column("barcode") %>%
  mutate(tmp.cluster = ifelse(
    harmony.cluster.0.4 %in% cluster1, "cluster1", "cluster2"
  )) %>%
  column_to_rownames("barcode")
tmp.metadata <- tmp.metadata[row.names(s.obj.subset@meta.data), ]
s.obj.subset <- AddMetaData(object = s.obj.subset, metadata = tmp.metadata$tmp.cluster, col.name = "tmp.cluster")


DimPlot(object = s.obj.subset, reduction = "harmony_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "tmp.cluster")

tmp.umapdf <- s.obj.subset@reductions$harmony_UMAP@cell.embeddings %>% 
  as.data.frame() %>% 
  rownames_to_column("barcode")
tmp.umapdf <- merge(tmp.umapdf, s.obj.subset@meta.data %>% rownames_to_column("barcode"), by.x = "barcode", by.y = "barcode")

remove.cells <- subset(
  tmp.umapdf, tmp.umapdf$tmp.cluster == "cluster2" & tmp.umapdf$harmonyUMAP_1 <= 5
)$barcode

s.obj.subset <- subset(s.obj.subset, cells = setdiff(colnames(s.obj.subset), remove.cells))

p <- DimPlot(object = s.obj.subset, reduction = "harmony_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "tmp.cluster")
ggsave(plot = p, filename = sprintf("%s_%s.rename.pdf", selected.celltype, selected.label), 
       path = file.path(path.to.save.output, selected.label, selected.celltype), 
       device = "pdf", 
       width = 14, 
       height = 10, 
       dpi = 300)

##### perform DGE analysis
if (file.exists(file.path(path.to.save.output, selected.label, selected.celltype, "diff_markers.xlsx")) == FALSE){
  s.obj.subset <- PrepSCTFindMarkers(s.obj.subset)
  diff.markers <- FindMarkers(object = s.obj.subset, ident.1 = "cluster1", ident.2 = "cluster2", group.by = "tmp.cluster", assay = "SCT", test.use = "wilcox")
  diff.markers <- diff.markers %>% as.data.frame() %>% rownames_to_column("Gene") %>%#
    subset(p_val_adj <= 0.05) %>%
    rowwise() %>%
    mutate(abs_avg_log2FC = abs(avg_log2FC)) %>%
    arrange(desc(avg_log2FC))
  
  writexl::write_xlsx(diff.markers, file.path(path.to.save.output, selected.label, selected.celltype, "diff_markers.xlsx"))  
}

##### plot some marker genes
plot.genes <- list(
  Astrocytes = c(
    "AQP4",
     "GFAP",
     "SLC1A2",
     "SLC1A3",
     "MLC1",
     "HEPACAM",
     "APOE",
     "GLUL",
     "KCNJ10"
    ),
  Oligodendrocytes = c(
    "PLP1",
    "OPALIN",
    "MOBP",
    "GPM6A",
    "MAG",
    "ERMN",
    "CLDN11",
    "UGT8",
    "NFASC"
  )
)

for (group in names(plot.genes)){
  print(sprintf("working on group %s", group))
  dir.create(file.path(path.to.save.output, "some_marker_genes", group), showWarnings = FALSE, recursive = TRUE)
  
  for (input.gene in plot.genes[[group]]){
    print(sprintf("working on input gene %s", input.gene))
    p.full <- FeaturePlot(object = s.obj.integrated, reduction = "harmony_UMAP", features = input.gene, label = TRUE)
    ggsave(plot = p.full, filename = sprintf("%s.full.pdf", input.gene), 
           path = file.path(path.to.save.output, "some_marker_genes", group), 
           device = "pdf", 
           width = 14, 
           height = 10, 
           dpi = 300)
    
    for (input.label in unique(s.obj.integrated$label)){
      print(sprintf("working on subset %s", input.label))
      tmp.subset.s.obj <- subset(s.obj.integrated, label == input.label)
      p.subset <- FeaturePlot(object = tmp.subset.s.obj, reduction = "harmony_UMAP", features = input.gene, label = TRUE)
      ggsave(plot = p.subset, filename = sprintf("%s.%s.pdf", input.gene, input.label), 
             path = file.path(path.to.save.output, "some_marker_genes", group), 
             device = "pdf", 
             width = 14, 
             height = 10, 
             dpi = 300)
    }
  }
}


