gc()
rm(list = ls())

path.to.main.src <- "/media/hieunguyen/HNSD01/src/UKK/src/BrainMet/scRNAseq"
path.to.rmd <- file.path(path.to.main.src, "01_preliminary_analysis.Rmd")

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))

outdir <- "/media/hieunguyen/HNSD_mini/data/outdir"
input.outdir <- "/media/hieunguyen/HNHD01/outdir"

project.type <- "Lung"

fetch_all_folders <- basename(Sys.glob(file.path(input.outdir, "Lung*")))

all.PROJECTS <- to_vec( for(item in fetch_all_folders) sprintf("%s_%s", str_split(item, "_SeuratV5_")[[1]][[1]], "SeuratV5")) %>% unique()
for (PROJECT in all.PROJECTS){
  all.datasets <- to_vec( for (item in fetch_all_folders) if(grepl(PROJECT, item) == TRUE) str_split(str_replace(item, sprintf("%s_", PROJECT), ""), "_")[[1]][[1]]) %>% unique()
  for (dataset.name in all.datasets){
    ##### adjust input params for some special cases
    if (dataset.name == "p019n"){
      dataset.name <- "p019n_p028n_p030n_p033n"
    } else if (dataset.name == "KU"){
      dataset.name <- "KU_loom"
    }
    
    all.samples <- to_vec (for (item in fetch_all_folders) if(grepl(sprintf("%s_%s", PROJECT, dataset.name), item) == TRUE ) str_replace(item, sprintf("%s_", PROJECT), ""))
    for (sample.id in all.samples){
      print(sprintf("PROJECT: %s, dataset: %s, sample: %s", PROJECT, dataset.name, sample.id))
      path.to.save.html <- file.path(outdir,"html_outputs", PROJECT, dataset.name)
      dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)
      save.html.filename <- sprintf("01_preliminary_analysis_%s_%s_%s.html", PROJECT, dataset.name, sample.id)
      if (file.exists(file.path(path.to.save.html, save.html.filename)) == FALSE){
        rmarkdown::render(input = path.to.rmd, 
                          output_dir = path.to.save.html, 
                          output_file = save.html.filename, 
                          params = list(outdir = outdir, 
                                        input.outdir = input.outdir,
                                        PROJECT = PROJECT,
                                        dataset.name = dataset.name,
                                        sample.id = sample.id,
                                        project.type = project.type))       
        ##### move output to HDD storage, release working space on SSD storage
        dir.create(file.path(input.outdir, sprintf("%s_%s", PROJECT, sample.id), sprintf("%s", sample.id), "data_analysis"))
        system(sprintf("mv  %s/* %s", 
                       file.path(outdir, PROJECT, dataset.name, sample.id),
                       file.path(input.outdir, sprintf("%s_%s", PROJECT, sample.id), sprintf("%s", sample.id), "data_analysis")))
        system(sprintf("rm -rf %s", file.path(outdir, PROJECT, dataset.name, sample.id)))
      } else {
        print(sprintf("File %s exists!", file.path(path.to.save.html, save.html.filename)))
      }
    }
  }
}

