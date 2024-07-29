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

outdir <- "/media/hieunguyen/HNSD01/outdir"
main.PROJECT <- "BrainMet_SeuratV5"

path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))
path.to.10.output <- file.path(path.to.main.output, "10_output")
dir.create(path.to.10.output, showWarnings = FALSE, recursive = TRUE)

maindir <- "/media/hieunguyen/HNHD01/outdir/BrainMet_SeuratV5_GSE193745_GSE193745/GSE193745_GSE193745/s3_output"
main.metadatadir <- "/media/hieunguyen/HNHD01/outdir/BrainMet_SeuratV5_GSE193745_GSE193745/GSE193745_GSE193745/s8a_output"
savedir <- "/media/hieunguyen/HNHD01/outdir"

s.obj <- readRDS(file.path(maindir, "BrainMet_SeuratV5_GSE193745_GSE193745.output.s3.rds"))
download.meta.data <- read.csv(file.path(main.metadatadir, "GSE193745_Mets2022.TIL.metadata.txt"), sep = "\t") %>%
  rownames_to_column("barcode")
print("finished reading input data to s.obj")

meta.data <- s.obj@meta.data %>% rownames_to_column("rowname") %>%
  mutate(barcode = str_replace(rowname, "GSE193745_GSE193745_", ""))

meta.data <- merge(meta.data, subset(download.meta.data, select = c(barcode, ID, CancerType)), by.x = "barcode", by.y = "barcode")

meta.data <- meta.data %>% rowwise() %>%
  mutate(name = ifelse(grepl("_", ID) == TRUE, str_split(ID, "_")[[1]][[1]], ID)) %>%
  column_to_rownames("rowname")

meta.data <- meta.data[row.names(s.obj@meta.data), ]
s.obj <- AddMetaData(object = s.obj, metadata = meta.data$CancerType, col.name = "PrimaryTumor")
s.obj <- AddMetaData(object = s.obj, metadata = meta.data$name, col.name = "name")
print("finished adding metadata sample name for s.obj")

sample.id <- "merge17samples"
dir.create(file.path(savedir, sprintf("BrainMet_SeuratV5_GSE193745_%s", sample.id), sprintf("GSE193745_%s", sample.id), "s3_output"), showWarnings = FALSE, recursive = TRUE)

print("saving file ...")
saveRDS(s.obj, 
        file.path(savedir, 
                  sprintf("BrainMet_SeuratV5_GSE193745_%s", sample.id), 
                  sprintf("GSE193745_%s", sample.id), 
                  "s3_output", 
                  sprintf("BrainMet_SeuratV5_GSE193745_%s.output.s3.rds", sample.id)))
print('finished')
# if (file.exists(file.path(maindir, "finished_splitting_GSE193745.csv")) == FALSE){
#   ##### NOTE
#   # split the big dataset into sample-wise dataset
#   # save to the same folder with the same structure, which can be loaded easily 
#   # by the merge object function in 02 and integrate in 03.
#   
#   # One of the sample, LB3935T, was excluded from the combined counts matrix and 
#   # subsequent integrated analysis, due to low cell number (<=100 cells).
#   # link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193745
#   #####
#   dataset.metadata <- data.frame()
#   for (sample.id in unique(s.obj$name)){
#     dir.create(file.path(savedir, sprintf("BrainMet_SeuratV5_GSE193745_%s", sample.id), sprintf("GSE193745_%s", sample.id), "s3_output"), showWarnings = FALSE, recursive = TRUE)
#     print(sprintf("saving %s ...", sample.id))
#     tmp.s.obj <- subset(s.obj, name == sample.id)
#     saveRDS(tmp.s.obj, file.path(savedir, sprintf("BrainMet_SeuratV5_GSE193745_%s", sample.id), sprintf("GSE193745_%s", sample.id), "s3_output", 
#                                  sprintf("BrainMet_SeuratV5_GSE193745_%s.output.s3.rds", sample.id)))
#     
#     tmp.dataset.metadata <- data.frame(PROJECT = "BrainMet_SeuratV5",
#                                        Dataset = "GSE193745",
#                                        Sample = sample.id,
#                                        Group = "BrainMet")
#     dataset.metadata <- rbind(dataset.metadata, tmp.dataset.metadata)
#   }
#   write.csv(dataset.metadata, file.path(maindir, "GSE193745_sample_metadata.csv"))
#   write.csv(data.frame(status = c("finished splitting data GSE193745")), file.path(maindir, "finished_splitting_GSE193745.csv"))
# }
