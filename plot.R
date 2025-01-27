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

path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))
path.to.03.output <- file.path(path.to.main.output, "03_output")
path.to.09.output <- file.path(path.to.main.output, "09_output")
path.to.10.output <- file.path(path.to.main.output, "10_output")
path.to.11.output <- file.path(path.to.main.output, "11_output")
path.to.12.output <- file.path(path.to.main.output, "12_output")

s.obj.integrated <- readRDS(file.path(path.to.12.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds"))

Idents(s.obj.integrated) <- "harmony.cluster.0.5"

s.obj.cluster21 <- subset(s.obj.integrated, harmony.cluster.0.5 == 21)

count.sample.in.cluster21 <- table(s.obj.cluster21$name) %>% as.data.frame()
writexl::write_xlsx(count.sample.in.cluster21, file.path(path.to.12.output, "num_cells_from_datasets_in_cluster_21.xlsx"))

GSE193745.metadata.path <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq/integrate_BrainMet_datasets/GSE193745_Mets2022.TIL.metadata.txt"
GSE193745.metadata <- read.csv(GSE193745.metadata.path, sep = "\t") %>% rownames_to_column("barcode") %>% subset(select = c(barcode, ID, CancerType))

s.obj.merge17 <- subset(s.obj.integrated, name == "merge17samples") 

merge17.metadata <- s.obj.merge17@meta.data %>% rownames_to_column("barcode") %>% as.data.frame() %>%
  rowwise() %>%
  mutate(barcode = str_replace(barcode, "GSE193745_GSE193745_", ""))

merge17.metadata <- merge(merge17.metadata, GSE193745.metadata, by.x = "barcode", by.y = "barcode")

merge17.metadata.cluster21 <- subset(merge17.metadata, merge17.metadata$harmony.cluster.0.5 == 21)

count.cancer.type <- table(merge17.metadata.cluster21$CancerType) %>% as.data.frame()
writexl::write_xlsx(count.cancer.type, file.path(path.to.12.output, "num_cells_cancer_type_cluster21_GSE193745.xlsx"))
count.samples <- table(merge17.metadata.cluster21$ID) %>% as.data.frame()
writexl::write_xlsx(count.samples, file.path(path.to.12.output, "num_cells_samples_cluster21_GSE193745.xlsx"))

##### Update 19.11.2024
gene.list <- c("CD3E", "CD4", "CD8A", "NCAM1", "ITGAX", "FCGR3B", "HLA-DRA")
DefaultAssay(s.obj.integrated) <- "SCT"
feature.plot <- FeaturePlot(object = s.obj.integrated, reduction = "harmony_UMAP", features = gene.list, label = TRUE, pt.size = 1) &
  scale_color_gradient(low = "lightgray", high = "#FF0000", na.value = "lightgray")
# scale_color_distiller(palette = "RdBu")
violin.plot <- VlnPlot(object = s.obj.integrated, features = gene.list, pt.size = 0, group.by = "harmony.cluster.0.5")
s.obj.integrated <- AddModuleScore(object = s.obj.integrated, features = list(gene.list), name = "check_gene_list_", ctrl = 50)

feature.plot.module <- FeaturePlot(object = s.obj.integrated, reduction = "harmony_UMAP", features = "check_gene_list_1", label = TRUE, pt.size = 1) &
  scale_color_gradient(low = "lightgray", high = "#FF0000", na.value = "lightgray")
violin.plot.module <- VlnPlot(object = s.obj.integrated, features = "check_gene_list_1", pt.size = 0, group.by = "harmony.cluster.0.5")

