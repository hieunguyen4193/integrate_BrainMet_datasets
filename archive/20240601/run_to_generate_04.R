gc()
rm(list = ls())

path.to.main.src <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq"
path.to.rmd <- file.path(path.to.main.src, "01_preliminary_analysis.Rmd")

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

outdir <- "/media/hieunguyen/HNSD01/outdir"
output.version <- "20240601"
path.to.rmd <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq/integrate_BrainMet_datasets/04_downstream_analysis.Rmd"

integrated.version <- "v0.2"
outdir <- "/media/hieunguyen/HNSD01/outdir"
main.PROJECT <- "BrainMet_SeuratV5"
code.version <- "20240601"
PROJECT <- "integrated_BrainMet_dataset"

path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")
path.to.06.output <- file.path(path.to.main.output, "06_output")

for (input.path.to.sobj in c(
  file.path(path.to.06.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)),
  file.path(path.to.03.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT))
)){
  path.to.save.html <- file.path(outdir, "BrainMet_SeuratV5", output.version, "html_outputs")
  dir.create(file.path(path.to.save.html), showWarnings = FALSE, recursive = TRUE)
  if (grepl("03_output", input.path.to.sobj) == TRUE){
    save.html.name <- "all_cells.html"
  } else if (grepl("06_output", input.path.to.sobj) == TRUE){
    save.html.name <- "immune_cells_CD45_only.html"
  }
  if (file.exists(file.path(path.to.save.html, save.html.name)) == FALSE){
    rmarkdown::render(path.to.rmd, 
                      params = list( path.to.sobj = input.path.to.sobj),
                      output_file = save.html.name,
                      output_dir = path.to.save.html)
  }
}
