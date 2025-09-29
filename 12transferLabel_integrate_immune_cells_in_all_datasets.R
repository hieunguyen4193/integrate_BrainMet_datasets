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

output.version <- "20250928"
# recently added, keep many versions. 
path.to.save.output <- file.path(path.to.main.output, "transfer_label_from_BrainMet", output.version)
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

# perform standard preprocessing on each object
ref.obj <- all.s.obj$BrainMet

predicted.metadata <- list()

for (input.label in c("Control Epilepsy", "Control Glioma")){
  if (file.exists(file.path(path.to.save.output, sprintf("%s.transferAnnotation.csv", input.label))) == FALSE){
    query.obj <- all.s.obj[[input.label]]
    DefaultAssay(ref.obj) <- "RNA"
    DefaultAssay(query.obj) <- "RNA"
    
    ref.obj <- NormalizeData(ref.obj)
    ref.obj <- FindVariableFeatures(ref.obj)
    ref.obj <- ScaleData(ref.obj)
    
    query.obj <- NormalizeData(query.obj)
    query.obj <- FindVariableFeatures(query.obj)
    query.obj <- ScaleData(query.obj)
    
    anchors <- FindTransferAnchors(reference = ref.obj, query = query.obj)
    
    predictions <- TransferData(
      anchorset = anchors,
      refdata = ref.obj$celltype
    )
    query.obj <- AddMetaData(object = query.obj, metadata = predictions, col.name = "predicted.id")
    
    DefaultAssay(query.obj) <- "SCT"
    p <- DimPlot(object = query.obj, 
                 label = TRUE, 
                 label.box = TRUE, 
                 group.by = "predicted.id", 
                 reduction = "harmony_UMAP")
    save.metadata <- query.obj@meta.data %>% rownames_to_column("barcode") %>% subset(select = c(barcode, predicted.id))
    write.csv(save.metadata, 
              file.path(path.to.save.output, sprintf("%s.transferAnnotation.csv", input.label)))
    predicted.metadata[[input.label]] <- save.metadata
  } else {
    print("reading in saved metadata")
    predicted.metadata[[input.label]] <- read.csv(file.path(path.to.save.output, sprintf("%s.transferAnnotation.csv", input.label)))
  }
} 