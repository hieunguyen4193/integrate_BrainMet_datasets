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
outdir <- "/media/hieunguyen/HNSD01/outdir"
main.PROJECT <- "BrainMet_SeuratV5"

sample.metadata <- all.integrated.metadata[[integrated.version]]

main.PROJECT <- "BrainMet_SeuratV5"
code.version <- "integrate_BrainMet_datasets"
PROJECT <- "integrated_BrainMet_dataset"

path.to.main.output <- file.path(outdir, main.PROJECT, code.version, sprintf("integrated_%s", integrated.version))

output.version <- '20240730'
path.to.save.html <- file.path(path.to.main.output, "html_outputs", output.version)
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

path.to.03.output <- file.path(path.to.main.output, "03_output")
path.to.09.output <- file.path(path.to.main.output, "09_output")
path.to.10.output <- file.path(path.to.main.output, "10_output")
path.to.11.output <- file.path(path.to.main.output, "11_output")
path.to.12.output <- file.path(path.to.main.output, "12_output")
path.to.13.output <- file.path(path.to.main.output, "13_output")

input.list <- list(
  v2_all_cells = list(
    path.to.s.obj = file.path(path.to.03.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.5),
  v2_immune_cells_0.5 = list(
    path.to.s.obj = file.path(path.to.main.output, "06_output_0.5", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.5
  ),
  v2_immune_cells_0.75 = list(
    path.to.s.obj = file.path(path.to.main.output, "06_output_0.75", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.75
  ),
  v2_immune_cells_1 = list(
    path.to.s.obj = file.path(path.to.main.output, "06_output_1", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 1
  ),
  v2_immune_cells_Control_Epilepsy_0.5 = list(
    path.to.s.obj = file.path(path.to.main.output, "07_output_Control Epilepsy_0.5", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.5
  ),
  v2_immune_cells_Control_Epilepsy_0.75 = list(
    path.to.s.obj = file.path(path.to.main.output, "07_output_Control Epilepsy_0.75", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.75
  ),
  v2_immune_cells_Control_Epilepsy_1 = list(
    path.to.s.obj = file.path(path.to.main.output, "07_output_Control Epilepsy_1", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 1
  ),
  v2_immune_cells_Control_Glioma_0.5 = list(
    path.to.s.obj = file.path(path.to.main.output, "07_output_Control Glioma_0.5", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.5
  ),
  v2_immune_cells_Control_Glioma_0.75 = list(
    path.to.s.obj = file.path(path.to.main.output, "07_output_Control Glioma_0.75", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.75
  ),
  v2_immune_cells_Control_Glioma_1 = list(
    path.to.s.obj = file.path(path.to.main.output, "07_output_Control Glioma_1", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 1
  ),
  v2_immune_cells_BrainMet_0.5 = list(
    path.to.s.obj = file.path(path.to.main.output, "07_output_BrainMet_0.5", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.5
  ),
  v2_immune_cells_BrainMet_0.75 = list(
    path.to.s.obj = file.path(path.to.main.output, "07_output_BrainMet_0.75", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.75
  ),
  v2_immune_cells_BrainMet_1 = list(
    path.to.s.obj = file.path(path.to.main.output, "07_output_BrainMet_1", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 1
  ),
  v2_filtered_immune_cells_BrainMet_0.5 = list(
    path.to.s.obj = file.path(path.to.main.output, "08_output_BrainMet_0.5", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.5
  ),
  v2_filtered_immune_cells_Control_Epilepsy_0.75 = list(
    path.to.s.obj = file.path(path.to.main.output, "08_output_Control Epilepsy_0.75", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.75
  ),
  v2_filtered_immune_cells_Control_Glioma_0.5 = list(
    path.to.s.obj = file.path(path.to.main.output, "08_output_Control Glioma_0.5", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.5
  ),
  v2_filtered_immune_cells_0.5 = list(
    path.to.s.obj = file.path(path.to.main.output, "09_output", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.5
  ),
  v2_GSE193745_integrated_0.5 = list(
    path.to.s.obj = file.path(path.to.main.output, "10_output", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.5
  ),
  GSE193745_immune_cells_integrated_0.5 = list(
    path.to.s.obj = file.path(path.to.main.output, "11_output", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.5
  ),
  v2_and_GSE193745_immune_cells_integrated_0.5 = list(
    path.to.s.obj = file.path(path.to.main.output, "12_output", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.5
  ),
  v2_and_GSE193745_all_cells_integrated_0.5 = list(
    path.to.s.obj = file.path(path.to.main.output, "13_output", "s8_output", "integrated_BrainMet_dataset.output.s8.rds"),
    cluster.resolution = 0.5
  )
)

for (input.case in names(input.list)){
  input.path.to.sobj <- input.list[[input.case]]$path.to.s.obj
  cluster.resolution <- input.list[[input.case]]$cluster.resolution
  outputdir <- dirname(input.path.to.sobj)
  save.html.name <- sprintf("%s.html", input.case)
  if (file.exists(file.path(path.to.save.html, save.html.name)) == FALSE){
    rmarkdown::render(path.to.rmd,
                      params = list(
                        path.to.sobj = input.path.to.sobj,
                        cluster.resolution = cluster.resolution,
                        outputdir = outputdir),
                      output_file = save.html.name,
                      output_dir = path.to.save.html)
  }
}
