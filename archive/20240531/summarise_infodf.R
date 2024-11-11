gc()
rm(list = ls())


all.infodf <- Sys.glob(file.path("/media/hieunguyen/HNHD01/outdir/*/*/data_analysis/01_output/infodf.xlsx"))
infodf <- data.frame()
for (file in all.infodf){
  infodf <- rbind(infodf, readxl::read_excel(file))
}

writexl::write_xlsx(infodf, file.path("/media/hieunguyen/HNSD_mini/data/outdir/html_outputs/info.xlsx"))

project.id <- "BrainMet_SeuratV5"

infodf <- subset(infodf, infodf$PROJECT == project.id)

table(infodf$Species)
sum(infodf$Num_cells)
