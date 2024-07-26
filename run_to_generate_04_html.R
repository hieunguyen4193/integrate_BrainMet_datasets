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
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")

#####----------------------------------------------------------------------#####
##### Running 06 and 07 data analysis
#####----------------------------------------------------------------------#####
# for (input.condition in unique(sample.metadata$Group)){
#   for (cluster.resolution in c(0.5, 0.75, 1)){
#     # path.to.06.output <- file.path(path.to.main.output, sprintf("06_output_%s", cluster.resolution))
#     # input.path.to.sobj <- file.path(path.to.06.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds")
#     path.to.07.output <- file.path(path.to.main.output, sprintf("07_output_%s_%s", input.condition, cluster.resolution))
#     input.path.to.sobj <- file.path(path.to.07.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds")
#     path.to.save.html <- file.path(outdir, "BrainMet_SeuratV5", output.version, "html_outputs")
#     dir.create(file.path(path.to.save.html), showWarnings = FALSE, recursive = TRUE)
#     if (grepl("03_output", input.path.to.sobj) == TRUE){
#       save.html.name <- "all_cells.html"
#     } else if (grepl("06_output", input.path.to.sobj) == TRUE){
#       save.html.name <- sprintf("immune_cells_CD45_only.clusterRes_%s.html", cluster.resolution)
#     } else if (grepl("07_output", input.path.to.sobj) == TRUE){
#       save.html.name <- sprintf("immune_cells_CD45_only_condition_%s.clusterRes_%s.html", input.condition, cluster.resolution)
#     }
#     if (file.exists(file.path(path.to.save.html, save.html.name)) == FALSE){
#       rmarkdown::render(path.to.rmd, 
#                         params = list( 
#                           path.to.sobj = input.path.to.sobj,
#                           cluster.resolution = cluster.resolution),
#                         output_file = save.html.name,
#                         output_dir = path.to.save.html)
#     }
#   }
# }

#####----------------------------------------------------------------------#####
##### running data analysis 04 template for 09 output data
#####----------------------------------------------------------------------#####
# path.to.09.output <- file.path(path.to.main.output, "09_output")
# 
# cluster.resolution <- 0.5
# input.path.to.sobj <- file.path(path.to.09.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds")
# 
# path.to.save.html <- file.path(outdir, "BrainMet_SeuratV5", output.version, "html_outputs")
# dir.create(file.path(path.to.save.html), showWarnings = FALSE, recursive = TRUE)
# 
# if (grepl("03_output", input.path.to.sobj) == TRUE){
#   save.html.name <- "all_cells.html"
# } else if (grepl("06_output", input.path.to.sobj) == TRUE){
#   save.html.name <- sprintf("immune_cells_CD45_only.clusterRes_%s.html", cluster.resolution)
# } else if (grepl("07_output", input.path.to.sobj) == TRUE){
#   save.html.name <- sprintf("immune_cells_CD45_only_condition_%s.clusterRes_%s.html", input.condition, cluster.resolution)
# } else if (grepl("09_output", input.path.to.sobj) == TRUE){
#   save.html.name <- sprintf("immune_cells_CD45_only.clusterRes_%s.filtered.html", cluster.resolution)
# }
# if (file.exists(file.path(path.to.save.html, save.html.name)) == FALSE){
#   rmarkdown::render(path.to.rmd,
#                     params = list(
#                       path.to.sobj = input.path.to.sobj,
#                       cluster.resolution = cluster.resolution),
#                     output_file = save.html.name,
#                     output_dir = path.to.save.html)
# }


##### running data analysis 04 template for 10 output data

path.to.10.output <- file.path(path.to.main.output, "09_output")

cluster.resolution <- 0.5
input.path.to.sobj <- file.path(path.to.10.output, "s8_output", "integrated_BrainMet_dataset.output.s8.rds")

path.to.save.html <- file.path(outdir, "BrainMet_SeuratV5", output.version, "html_outputs")
dir.create(file.path(path.to.save.html), showWarnings = FALSE, recursive = TRUE)

if (grepl("03_output", input.path.to.sobj) == TRUE){
  save.html.name <- "all_cells.html"
} else if (grepl("06_output", input.path.to.sobj) == TRUE){
  save.html.name <- sprintf("immune_cells_CD45_only.clusterRes_%s.html", cluster.resolution)
} else if (grepl("07_output", input.path.to.sobj) == TRUE){
  save.html.name <- sprintf("immune_cells_CD45_only_condition_%s.clusterRes_%s.html", input.condition, cluster.resolution)
} else if (grepl("09_output", input.path.to.sobj) == TRUE){
  save.html.name <- sprintf("immune_cells_CD45_only.clusterRes_%s.filtered.html", cluster.resolution)
} else if (grepl("10_output", input.path.to.sobj) == TRUE){
  save.html.name <- sprintf("immune_cells_CD45_only.clusterRes_%s.filtered.html", cluster.resolution)
}
if (file.exists(file.path(path.to.save.html, save.html.name)) == FALSE){
  rmarkdown::render(path.to.rmd,
                    params = list(
                      path.to.sobj = input.path.to.sobj,
                      cluster.resolution = cluster.resolution),
                    output_file = save.html.name,
                    output_dir = path.to.save.html)
}


