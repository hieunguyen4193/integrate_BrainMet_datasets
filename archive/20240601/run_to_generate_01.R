gc()
rm(list = ls())

path.to.main.src <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq"
path.to.rmd <- file.path(path.to.main.src, "01_preliminary_analysis.Rmd")

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

outdir <- "/media/hieunguyen/HNSD_mini/data/outdir"
input.outdir <- "/media/hieunguyen/HNHD01/outdir"

sample.metadata <- read.csv("/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq/20240601/samples_to_integrated_20240513.csv")
output.version <- "20240612"
path.to.rmd <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq/20240601/01_preliminary_analysis.Rmd"

for (i in seq(1, nrow(sample.metadata))){
  PROJECT <- sample.metadata[i, ]$PROJECT
  dataset.name <- sample.metadata[i, ]$Dataset
  sample.id <- sample.metadata[i, ]$Sample
  input.outdir <- "/media/hieunguyen/HNHD01/outdir"
  project.type <- "Brain"
  
  path.to.save.html <- file.path(outdir, "BrainMet_SeuratV5", output.version, "html_outputs")
  dir.create(file.path(path.to.save.html), showWarnings = FALSE, recursive = TRUE)
  save.html.name <- sprintf("%s_%s.html", dataset.name, sample.id)
  if (file.exists(file.path(path.to.save.html, save.html.name)) == FALSE){
    rmarkdown::render(path.to.rmd, params = list( PROJECT = PROJECT,
                                                  outdir = outdir,
                                                  dataset.name = dataset.name,
                                                  sample.id= sample.id,
                                                  input.outdir= input.outdir,
                                                  project.type = project.type),
                      output_file = save.html.name,
                      output_dir = path.to.save.html)
  }
}
