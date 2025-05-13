#####----------------------------------------------------------------------#####
##### 01
#####----------------------------------------------------------------------#####
# gc()
# rm(list = ls())

convert_Seurat_to_anndata <- function(path.to.s.obj, outputdir, umap.name, pca.name, cluster.name){
  s.obj <- readRDS(path.to.s.obj)
  PROJECT <- str_replace(basename(path.to.s.obj), ".rds", "")
  
  outputdir <- file.path(outputdir, PROJECT)
  dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
  
  s.obj$barcode <- colnames(s.obj)
  
  s.obj$UMAP_1 <- s.obj@reductions[[umap.name]]@cell.embeddings[,1]
  s.obj$UMAP_2 <- s.obj@reductions[[umap.name]]@cell.embeddings[,2]
  
  write.csv(s.obj@reductions[[pca.name]]@cell.embeddings, 
            file=file.path(outputdir, 
                           sprintf('pca_%s.csv', sprintf("%s_%s", PROJECT, cluster.name))), 
            quote=F, 
            row.names=F)
  
  write.csv(s.obj@meta.data, file=file.path(outputdir, sprintf('metadata_%s.csv', sprintf("%s_%s", PROJECT, cluster.name))), quote=F, row.names=F)
  
  # write expression counts matrix
  library(Matrix)
  counts_matrix <- GetAssayData(s.obj, assay='SCT', slot='data')
  writeMM(counts_matrix, file=file.path(outputdir, sprintf('counts_%s.mtx', sprintf("%s_%s", PROJECT, cluster.name))))
  
  # write gene names
  write.table( data.frame('gene'=rownames(counts_matrix)),file=file.path(outputdir, 
                                                                         sprintf('gene_names_%s.csv', 
                                                                                 sprintf("%s_%s", PROJECT, cluster.name))), quote=F,row.names=F,col.names=F)
  
}
